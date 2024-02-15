# RRR
# look at how many bins per dRep ANI 
# and how many bins from same sample combined per ANI
# to choose an ANI granularity to use.

library(ggplot2)
library(patchwork)
library(gridExtra)

drep.folder <- "data/2023-02-24_dRep_with_range_ANIs/output_processed"

output.folder <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/summaries"
plot.folder <- "figures/2023-02-27_dRep_ANI_comparisons"

ani.range <- seq(from = .90, to = .99, by = .01)

# ---- functions ----

make.empty.list.structure <- function(ListNames){
  empty.list <- list(NULL)
  for (e in 1:length(ListNames)){
    empty.list[[e]] <- 0
    names(empty.list)[e] <- ListNames[e]
  }
  return(empty.list)
}

get.clust.sum <- function(drep.list){
  clust.sum <- data.frame("ANI" = ani.range, "num.clusts" = NA, "num.combined" = NA)
  
  for (ani in ani.range){
    clust.sum$num.clusts[clust.sum$ANI == ani] <- sum(drep.list[[paste0("ANI_",ani)]]$winner)
  }
  
  for (ani in ani.range){
    my.d <- drep.list[[paste0("ANI_",ani)]]
    my.d <- aggregate(x = my.d[ ,paste0("drep.cluster.",ani,"ANI")], by = list("sample" = my.d$tymeflies.name, "cluster" = my.d[ ,paste0("drep.cluster.",ani,"ANI")]), FUN = length)
    clust.sum$num.combined[clust.sum$ANI == ani] <- sum(my.d$x > 1)
  }
  
  return(clust.sum)
}

subset.a.taxon <- function(tax.level, taxon, drep.list){
  for (ani in ani.range){
    my.d <- drep.list[[paste0("ANI_",ani)]]
    index <- which(my.d[ ,tax.level] == taxon)
    my.d <- my.d[index, ]
    drep.list[[paste0("ANI_",ani)]] <- my.d
  }
  return(drep.list)
}

plot.num.clusts <- function(my.sum, zoom = F){
  if(zoom){my.sum <- my.sum[1:9, ]}
  p <- ggplot(data = my.sum, aes(x = ANI, y = num.clusts))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    geom_point(color = "magenta4")+
    geom_line(color = "magenta4")+
    scale_x_continuous(name = element_blank(), breaks = ani.range)+
    scale_y_continuous(name = "Representative\nGenomes")
  return(p)
}

plot.num.combined <- function(my.sum, zoom = F){
  if(zoom){my.sum <- my.sum[1:9, ]}
  p <- ggplot(data = my.sum, aes(x = ANI, y = num.combined))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    geom_point(color = "cyan4")+
    geom_line(color = "cyan4")+
    scale_x_continuous(name = "dRep ANI", breaks = ani.range)+
    scale_y_continuous(name = "Same-sample\nBins Combined")+
    geom_hline(yintercept = 0)
  return(p)
}

plot.summary.table <- function(my.sum){
  
}

# ---- import data tables ----

drep.results <- make.empty.list.structure(ListNames = paste0("ANI_",ani.range))

drep.files <- list.files(path = drep.folder, pattern = ".rds", full.names = T)

for (ani in ani.range){
  drep.file <- grep(pattern = paste0(ani,".rds"), x = drep.files, value = T)
  drep.results[[paste0("ANI_",ani)]] <- readRDS(drep.file)
}

# ---- overall look ----

overall.summary <- get.clust.sum(drep.list = drep.results)
write.table(x = overall.summary, file = file.path(output.folder,"overall.tsv"), sep = "\t", quote = F, row.names = F)

p1 <- plot.num.clusts(my.sum = overall.summary)
p2 <- plot.num.combined(my.sum = overall.summary)
my.design <- "
AB
CB
"
pdf(file = file.path(plot.folder,"All_Bins.pdf"), width = 9, height = 4)
p1 + tableGrob(overall.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
  plot_layout(design = my.design, widths = c(2,1)) + 
  plot_annotation(title = "All Bins", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()

p1 <- plot.num.clusts(my.sum = overall.summary, zoom = T)
p2 <- plot.num.combined(my.sum = overall.summary, zoom = T)
pdf(file = file.path(plot.folder,"zoom_All_Bins.pdf"), width = 9, height = 4)
p1 + tableGrob(overall.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
  plot_layout(design = my.design, widths = c(2,1)) + 
  plot_annotation(title = "All Bins", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()
# ---- by-phylum looks ----

phyla <- unique(drep.results[[10]]$phylum)
phyla <- phyla[!is.na(phyla)]

for (phylum in phyla){
  phylum.bins <- subset.a.taxon(tax.level = "phylum", taxon = phylum, drep.list = drep.results)
  phylum.summary <- get.clust.sum(drep.list = phylum.bins)
  write.table(x = phylum.summary, file = file.path(output.folder,paste0(phylum,".tsv")), sep = "\t", quote = F, row.names = F)
  
  p1 <- plot.num.clusts(my.sum = phylum.summary)
  p2 <- plot.num.combined(my.sum = phylum.summary)
  pdf(file = file.path(plot.folder, paste0(phylum,".pdf")), width = 9, height = 4)
  p <- p1 + tableGrob(phylum.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
    plot_layout(design = my.design, widths = c(2,1)) + 
    plot_annotation(title = phylum, 
                    theme = theme(plot.title = element_text(hjust = .5)))
  print(p)
  dev.off()
  
  p1 <- plot.num.clusts(my.sum = phylum.summary, zoom = T)
  p2 <- plot.num.combined(my.sum = phylum.summary, zoom = T)
  pdf(file = file.path(plot.folder, paste0("zoom_",phylum,".pdf")), width = 9, height = 4)
  p <- p1 + tableGrob(phylum.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
    plot_layout(design = my.design, widths = c(2,1)) + 
    plot_annotation(title = phylum, 
                    theme = theme(plot.title = element_text(hjust = .5)))
  print(p)
  dev.off()
}

# ---- nanopelagicales looks ----

unique(drep.results[[10]]$genus[drep.results[[10]]$family == "f__Nanopelagicaceae"])

nano.bins <- subset.a.taxon(tax.level = "family", taxon = "f__Nanopelagicaceae", drep.list = drep.results)
nano.summary <- get.clust.sum(drep.list = nano.bins)
write.table(x = nano.summary, file = file.path(output.folder,"f__Nanopelagicaceae.tsv"), sep = "\t", quote = F, row.names = F)

b.bins <- subset.a.taxon(tax.level = "genus", taxon = "g__Nanopelagicus", drep.list = drep.results)
b.summary <- get.clust.sum(drep.list = b.bins)
write.table(x = b.summary, file = file.path(output.folder,"g__Nanopelagicus.tsv"), sep = "\t", quote = F, row.names = F)

a.bins <- subset.a.taxon(tax.level = "genus", taxon = "g__Planktophila", drep.list = drep.results)
a.summary <- get.clust.sum(drep.list = a.bins)
write.table(x = a.summary, file = file.path(output.folder,"g__Planktophila.tsv"), sep = "\t", quote = F, row.names = F)


p1 <- plot.num.clusts(my.sum = nano.summary)
p2 <- plot.num.combined(my.sum = nano.summary)
pdf(file = file.path(plot.folder,"f__Nanopelagicaceae.pdf"), width = 9, height = 4)
p1 + tableGrob(nano.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
  plot_layout(design = my.design, widths = c(2,1)) + 
  plot_annotation(title = "f__Nanopelagicaceae", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()

p1 <- plot.num.clusts(my.sum = nano.summary, zoom = T)
p2 <- plot.num.combined(my.sum = nano.summary, zoom = T)
pdf(file = file.path(plot.folder,"zoom_f__Nanopelagicaceae.pdf"), width = 9, height = 4)
p1 + tableGrob(nano.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
  plot_layout(design = my.design, widths = c(2,1)) + 
  plot_annotation(title = "f__Nanopelagicaceae", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()


p1 <- plot.num.clusts(my.sum = b.summary)
p2 <- plot.num.combined(my.sum = b.summary)
pdf(file = file.path(plot.folder,"g__Nanopelagicus.pdf"), width = 9, height = 4)
p1 + tableGrob(b.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
  plot_layout(design = my.design, widths = c(2,1)) + 
  plot_annotation(title = "g__Nanopelagicus", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()

p1 <- plot.num.clusts(my.sum = b.summary, zoom = T)
p2 <- plot.num.combined(my.sum = b.summary, zoom = T)
pdf(file = file.path(plot.folder,"zoom_g__Nanopelagicus.pdf"), width = 9, height = 4)
p1 + tableGrob(b.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
  plot_layout(design = my.design, widths = c(2,1)) + 
  plot_annotation(title = "g__Nanopelagicus", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()


p1 <- plot.num.clusts(my.sum = a.summary)
p2 <- plot.num.combined(my.sum = a.summary)
pdf(file = file.path(plot.folder,"g__Planktophila.pdf"), width = 9, height = 4)
p1 + tableGrob(a.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
  plot_layout(design = my.design, widths = c(2,1)) + 
  plot_annotation(title = "g__Planktophila", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()

p1 <- plot.num.clusts(my.sum = a.summary, zoom = T)
p2 <- plot.num.combined(my.sum = a.summary, zoom = T)
pdf(file = file.path(plot.folder,"zoom_g__Planktophila.pdf"), width = 9, height = 4)
p1 + tableGrob(a.summary, rows = NULL, cols = c("ANI","Number\nGenomes","Number\nCombined")) + p2 + 
  plot_layout(design = my.design, widths = c(2,1)) + 
  plot_annotation(title = "g__Planktophila", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()