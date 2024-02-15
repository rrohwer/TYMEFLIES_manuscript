# RRR
# get at whether same strains persist or new ones dominate
# from 2023-07-21 plot, but now remove low coverage samples. 

library(data.table)
library(vegan)
library(ggplot2)
library(lubridate)
library(patchwork)
library(viridisLite)

# local path testing
per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-21_3300043464_group3_bin69_SNVs.tsv.gz" # B
per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2016-07-20_3300033996_group7_bin32_SNVs.tsv.gz" # C
per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-04_3300044729_group3_bin142_SNVs.tsv.gz" # A

selection.summary.file <- "data/2023-10-26_selection_per_sample_summaries/example_data/ME2011-09-21_3300043464_group3_bin69_selection_summary.tsv.gz" # B
selection.summary.file <- "data/2023-10-26_selection_per_sample_summaries/example_data/ME2016-07-20_3300033996_group7_bin32_selection_summary.tsv.gz" # C
selection.summary.file <- "data/2023-10-26_selection_per_sample_summaries/example_data/ME2011-09-04_3300044729_group3_bin142_selection_summary.tsv.gz" # A

genome.info <- fread(input = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
genome.info[genome ==  "ME2011-09-21_3300043464_group3_bin69",above.breadth.cutoff := (breadth / breadth_expected) >= .5, ] # B
genome.info[genome ==  "ME2016-07-20_3300033996_group7_bin32",above.breadth.cutoff := (breadth / breadth_expected) >= .5, ] # C
genome.info[genome ==  "ME2011-09-04_3300044729_group3_bin142",above.breadth.cutoff := (breadth / breadth_expected) >= .5, ] # A

num.threads <- 1

# ---- 



# selection.summary <- fread(file = selection.summary.file, nThread = num.threads, colClasses = c("character","character","character","character","numeric","numeric","numeric"))

x <- "character"
names(x) <- "date"
snv.instrain <- fread(file = per.genome.snv.file, nThread = num.threads, colClasses = x)
my.genome <- snv.instrain[1, genome]

snv.instrain

# snv.instrain <- merge(x = snv.instrain, y = selection.summary[ ,c("sample","Npos","Nneg")], by = "sample")

genome.info <- genome.info[above.breadth.cutoff == TRUE & coverage_median > 10]

snv.instrain <- snv.instrain[sample %in% genome.info$sample]

snvs <- dcast(data = snv.instrain, formula = scaffold + position + gene + mutation_type ~ sample, value.var = "ref_freq", fill = 1)

snvs[1:10,1:7]
snvs$mutation_type[1:10] # note empty string when not on a gene!

snvs[ ,snv.num := paste0("snv",1:nrow(snvs))]
snv.key <- snvs[ ,.(scaffold, position, gene, mutation_type, snv.num)]
snv.key[ ,.(position = max(position)), by = scaffold] # curious how long the scaffolds are

snv.mat <- as.matrix(snvs[ ,-c(1:4)], rownames = "snv.num")
snv.mat <- t(snv.mat)

snv.dist.euc <- vegdist(x = snv.mat, method = "euclidean")
snv.nmds.euc <- metaMDS(comm = snv.dist.euc, trymax = 100)
stressplot(snv.nmds.euc)

# dates.key <- unique(snv.instrain[ ,.(sample, date, year, yday, season, invasion, Npos, Nneg)])
dates.key <- unique(snv.instrain[ ,.(sample, date, year, yday, season, invasion)])
# Try PCA instead

snv.pca <- prcomp(x = snv.mat, rank. = 2)
snv.pca$x[1:5, ]

my.pca <- data.table("sample" = rownames(snv.pca$x),"x" = snv.pca$x[ ,1], "y" = snv.pca$x[ ,2])
my.pca <- merge(x = my.pca, y = dates.key, by = "sample")
my.pca[ ,`:=`(season = factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")),
               invasion = factor(invasion, levels = c("none","spiny","zebra")),
               date = parse_date_time(x = date, orders = "ymd"))]
my.pca <- my.pca[order(date)]

# ---- make plots of long-term change ----

# euclidean:
my.nmds <- data.table("sample" = rownames(snv.nmds.euc$points),"x" = snv.nmds.euc$points[ ,1], "y" = snv.nmds.euc$points[ ,2])

my.nmds <- merge(x = my.nmds, y = dates.key, by = "sample")
my.nmds[ ,`:=`(season = factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")),
               invasion = factor(invasion, levels = c("none","spiny","zebra")),
               date = parse_date_time(x = date, orders = "ymd"))]
my.nmds <- my.nmds[order(date)]

# ---- for quick re-plotting ----


# under selection color scale

p.pos <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = Npos), alpha = .7, size = 3)+
  scale_color_viridis_c(option = "A")+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

p.pos

ggplot(data = my.pca, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = Npos), alpha = .7, size = 3)+
  scale_color_viridis_c(option = "A")+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

# annual color scale ----
p.drift <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = year), alpha = .7, size = 3)+
  scale_color_viridis_c(option = "D")+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

pdf(file = paste0("figures/2023-07-21_NMDS_plots_of_SNVs/plot_", my.genome,"_drift.pdf"), width = 5, height = 4)
p.drift
dev.off()

# add axis limits with coord_fixed to ensure exactly square for fair representation

# annual color scale with lines
p.drift.lines <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = year), alpha = .7, size = 3)+
  geom_path(aes(color = year), alpha = .3)+
  scale_color_viridis_c(option = "D")+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

pdf(file = paste0("figures/2023-07-21_NMDS_plots_of_SNVs/plot_", my.genome,"_drift_lines.pdf"), width = 5, height = 4)
p.drift.lines
dev.off()

# make years not continuous values
my.nmds[ ,year := factor(year)]  
my.colors <- c(rep("thistle3", 11), "purple2","red","green3","blue", rep("lightblue", 5))
names(my.colors) <- 2000:2019

# highlight change-years
p.2012 <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = year), alpha = .7, size = 2)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

pdf(file = paste0("figures/2023-07-21_NMDS_plots_of_SNVs/plot_", my.genome,"_2012.pdf"), width = 6.5, height = 6)
p.2012
dev.off()

# add lines connecting points
p.2012.lines <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = year), alpha = .7, size = 2)+
  geom_path(aes(color = year), alpha = .3)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

pdf(file = paste0("figures/2023-07-21_NMDS_plots_of_SNVs/plot_", my.genome,"_2012_lines.pdf"), width = 6.5, height = 6)
p.2012.lines
dev.off()

# color by season instead
my.colors <- c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")

# season colors
p.seasons <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season), alpha = .7)+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

pdf(file = paste0("figures/2023-07-21_NMDS_plots_of_SNVs/plot_", my.genome,"_seasons.pdf"), width = 5, height = 4)
p.seasons
dev.off()

# oof too much, color single-seasons:
# ICE
my.colors <- c("snow3","wheat","wheat","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(4,1,1,1,1,1)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.ice <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Ice-On")

# SPRING
my.colors <- c("wheat","tan4","wheat","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1,4,1,1,1,1)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.spring <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Spring")

# CLEAR
my.colors <- c("wheat","wheat","cornflowerblue","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1,1,4,1,1,1)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.clear <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Clearwater")

# EARLY
my.colors <- c("wheat","wheat","wheat","chartreuse4","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1,1,1,4,1,1)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.early <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Early Summer")

# LATE
my.colors <- c("wheat","wheat","wheat","wheat","purple","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1,1,1,1,4,1)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.late <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Late Summer")

# FALL
my.colors <- c("wheat","wheat","wheat","wheat","wheat","hotpink2")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1,1,1,1,1,4)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.fall <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Fall")

p.seasons.separate <- p.ice + p.spring + p.clear + p.early + p.late + p.fall + plot_layout(guides = "collect") +
  plot_annotation(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

pdf(file = paste0("figures/2023-07-21_NMDS_plots_of_SNVs/plot_", my.genome,"_seasons_separate.pdf"), width = 10, height = 8)
p.seasons.separate
dev.off()

# remove seasons from non-2012 years 
my.colors <- c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.nmds[ ,size := "big"]
my.nmds[year != 2012, `:=`(season = NA,
                           size = "small")]
my.nmds[ , size := factor(size, levels = c("small","big"))]

# but only color 2012 by season?
p.2012.seasons <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season, size = size), alpha = .7, show.legend = F)+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = c(1,4))+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)

pdf(file = paste0("figures/2023-07-21_NMDS_plots_of_SNVs/plot_", my.genome,"_2012_seasons.pdf"), width = 6.5, height = 6)
p.2012.seasons
dev.off()

my.nmds[ ,my.lab := date]
my.nmds[size == "small" ,my.lab := NA]
# get dates
ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = season, size = size), alpha = .7, show.legend = F)+
  geom_label(aes(label = my.lab))+
  geom_path(alpha = .1)+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = c(1,4))+
  labs(title = "Sample ordination based on SNV reference frequency", subtitle = my.genome)


# ---- end ----

