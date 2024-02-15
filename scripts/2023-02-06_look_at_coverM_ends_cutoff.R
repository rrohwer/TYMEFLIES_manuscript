# RRR

# ---- set-up ----

library(lubridate)
library(magrittr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

folder.path <- "data/2023-02-01_coverM_files/comparison_output/coverm_output/"

tax <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")

plot.folder.abs <- "figures/2023-02-01_coverM_end_cutoff_comparisons_absolute/"
plot.folder.perc <- "figures/2023-02-01_coverM_end_cutoff_comparisons_perc/"

# ---- functions ----


# ---- import and format dense ----

# my.files <- list.files(path = folder.path, pattern = "coverM_mean")
# 
# otu.list <- list(NULL)
# cutoffs <- NULL
# for (f in 1:length(my.files)){
#   otu <- read.table(file = file.path(folder.path,my.files[f]), header = T, sep = "\t")
#   row.names(otu) <- otu$Genome
#   otu <- otu[ ,-1]
#   otu <- as.matrix(otu)
#   cutoff <- sub(x = my.files[f], pattern = "coverM_mean_exclude_", replacement = "")
#   cutoff <- sub(x = cutoff, pattern = ".txt", replacement = "")
#   cutoffs[f] <- as.numeric(cutoff)
#   otu.list[[f]] <- otu
#   names(otu.list)[f] <- cutoff
# }
# 
# cutoffs
# 
# colnames(tax)
# 
# my.names <- data.frame("genome" = row.names(otu.list$`0`), "row.order" = 1:nrow(otu.list$`0`))
# 
# my.names <- merge(x = my.names, y = tax, by.x = "genome", by.y = "bin.full.name", all.x = TRUE, all.y = FALSE)
# my.names <- my.names[order(my.names$row.order), ]
# 
# sample.dates <- colnames(otu.list$`0`) %>%
#   sub(pattern = "^ME", replacement = "", x = .) %>%
#   sub(pattern = "_.*$", replacement = "", x = .) %>%
#   substr(x = ., start = 1, stop = 10) %>%
#   parse_date_time(orders = "ymd", tz = "Etc/GMT-5")

# ---- import and format long ----

my.files <- list.files(path = folder.path, pattern = "coverM_mean")

long.list <- list(NULL)
for (f in 1:length(my.files)){
  otu <- read.table(file = file.path(folder.path,my.files[f]), header = T, sep = "\t")
  otu <- pivot_longer(data = otu, cols = 2:ncol(otu), names_to = "sample.date", values_to = "mean")
  
  cutoff <- sub(x = my.files[f], pattern = "coverM_mean_exclude_", replacement = "") %>%
    sub(x = ., pattern = ".txt", replacement = "") %>%
    as.numeric()
  otu$cutoff <- cutoff
  
  otu$sample.date <- otu$sample.date %>%
    sub(pattern = "^ME", replacement = "", x = .) %>%
    sub(pattern = "_.*$", replacement = "", x = .) %>%
    substr(x = ., start = 1, stop = 10) %>%
    parse_date_time(orders = "ymd", tz = "Etc/GMT-5")
  
  long.list[[f]] <- otu
  names(long.list)[f] <- cutoff
}

long.table <- long.list[[1]]
for (e in 2:length(long.list)){
  long.table <- rbind(long.table, long.list[[e]])
}

rm(long.list)
rm(otu)
head(long.table)

colnames(tax)
colnames(long.table)

long.table <- merge(x = long.table, y = tax, by.x = "Genome", by.y = "bin.full.name", all.x = TRUE, all.y = FALSE)
rm(tax)

long.table <- long.table[order(long.table$Genome, long.table$cutoff, long.table$sample.date), ]
head(long.table)

# make it perc max coverage per day ----

take.perc <- function(x){
  x <- x / max(x)
  return(x)
}

max.daily.vals <- aggregate(x = long.table$mean, by = list(long.table$Genome, long.table$sample.date), FUN = max)
colnames(max.daily.vals) <- c("Genome","sample.date","max")

long.table <- merge(x = long.table, y = max.daily.vals, by = c("Genome","sample.date"), all = TRUE)

long.table$mean.perc <- long.table$mean / long.table$max * 100

overall.average <- aggregate(x = long.table[ ,c("mean","mean.perc")], by = list(long.table$cutoff), FUN = mean, na.rm = T)

phylum.average <- aggregate(x = long.table[ ,c("mean","mean.perc")], 
                            by = list("phylum" = long.table$phylum, "cutoff" = long.table$cutoff), FUN = mean, na.rm = T)

# ---- look at single genome single day first ----

# choose.genome <- 1
# choose.day <- 1
# 
# my.data <- data.frame("cutoff" = cutoffs, "mean" = NA)
# for (r in 1:nrow(my.data)){
#   my.data$mean[r] <- otu.list[[as.character(cutoffs[r])]][choose.genome,choose.day]
# }
# my.data <- my.data[order(my.data$cutoff), ]
# 
# ggplot(data = my.data, aes(x = cutoff, y = mean))+
#   geom_line()+
#   geom_point()

# ---- look at single genome all days overlaid ----

# # put it in a "long" format... kind of clunky
# 
# cutoffs <- cutoffs[order(cutoffs)]
# my.data <- data.frame("cutoff" = rep(rep(cutoffs, times = ncol(otu.list$`0`)), times = nrow(otu.list$`0`)), 
#                       "sample.date" = rep(colnames(otu.list$`0`), each = length(cutoffs), times = nrow(otu.list$`0`)), 
#                       "genome" = rep(row.names(otu.list$`0`), each = length(cutoffs) * ncol(otu.list$`0`)), 
#                       "mean" = NA)
# 
# for (r in 1:nrow(my.data)){
#   my.data$mean[r] <- otu.list[[as.character(my.data$cutoff[r])]][my.data$genome[r], my.data$sample.date[r]]
# }
# my.data <- my.data[order(my.data$cutoff), ]

# too slow, use tidyr to start with a long format:
# ----

my.genomes <- unique(long.table$Genome)
for (chosen.genome in my.genomes){
  index <- which(long.table$Genome == chosen.genome)
  genome.subset <- long.table[index, ]
  
  plot.tax <- paste(genome.subset$phylum[1],genome.subset$class[1],genome.subset$order[1],genome.subset$family[1],genome.subset$genus[1],genome.subset$species[1])
  
  p <- ggplot(data = genome.subset, aes(x = cutoff, y = mean, group = sample.date))+
    geom_line()+
    geom_point()+
    labs(title = chosen.genome, subtitle = plot.tax)

  png(filename = file.path(plot.folder.abs,paste(plot.tax, chosen.genome,".png")), width = 6, height = 3, units = "in", res = 70)
  print(p)
  dev.off()
  
  p.perc <-  ggplot(data = genome.subset, aes(x = cutoff, y = mean.perc, group = sample.date))+
    geom_line()+
    geom_point()+
    labs(title = chosen.genome, subtitle = plot.tax)
  
  png(filename = file.path(plot.folder.perc,paste(plot.tax, chosen.genome,".png")), width = 6, height = 3, units = "in", res = 70)
  print(p.perc)
  dev.off()
}

# ---- by-phylum summary ----


# ---- overall summary ----

overall.average

p.mean <- ggplot(data = overall.average, aes(x = Group.1, y = mean))+
  geom_line()+
  theme_bw()+
  geom_point()+
  geom_label(label = overall.average$Group.1, size = 3, label.padding = unit(0,"cm"), label.size = NA, 
             nudge_y = c(0,0,0,.001,.006,.006,.002,.002,.002,.002,.003,.005,.005,.006),
             nudge_x = c(-150,-150,-150,-175,-150,0,200,200,200,200,200,100,100,-25)) +
  labs(title = "Overall Average of all dRep Representative Bins")+
  scale_y_continuous(name = "Coverage (--methods mean)")+
  scale_x_continuous(name = "Length excluded (--contig-end-exclusion)")
  
p.perc <- ggplot(data = overall.average, aes(x = Group.1, y = mean.perc))+
  geom_line()+
  theme_bw()+
  geom_point()+
  geom_label(label = overall.average$Group.1, size = 2, label.padding = unit(0,"cm"), label.size = NA, 
             nudge_y = c(0,0,0,.25,.5,.75,.25,.25,.25,.25,.25,.25,.25,.25),
             nudge_x = c(-150,-150,-150,-175,-100,0,150,150,150,175,200,250,275,250)) +
  labs(title = "Overall Average of all dRep Representative Bins")+
  scale_y_continuous(name = "Percent Max Coverage (--methods mean)")+
  scale_x_continuous(name = "Length excluded (--contig-end-exclusion)")

png(filename = file.path(plot.folder.abs,"1-summary_overall_mean.png"), width = 5, height = 4, units = "in", res = 300)
p.mean
dev.off()

png(filename = file.path(plot.folder.perc,"1-summary_overall_mean-perc.png"), width = 5, height = 4, units = "in", res = 300)
p.perc
dev.off()
