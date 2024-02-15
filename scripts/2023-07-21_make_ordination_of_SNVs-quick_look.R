# RRR
# get at whether same strains persist or new ones dominate

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

num.threads <- 1

# ---- 

x <- "character"
names(x) <- "date"
snv.instrain <- fread(file = per.genome.snv.file, nThread = num.threads, colClasses = x)
my.genome <- snv.instrain[1, genome]

snv.instrain

snvs <- snv.instrain[ , .(scaffold, position, gene, mutation_type, ref_freq, sample) ]

snvs <- dcast(data = snvs, formula = scaffold + position + gene + mutation_type ~ sample, value.var = "ref_freq", fill = 1)

snvs[1:10,1:5]
snvs$mutation_type[1:10] # note empty string when not on a gene!

snvs[ ,snv.num := paste0("snv",1:nrow(snvs))]
snv.key <- snvs[ ,.(scaffold, position, gene, mutation_type, snv.num)]
snv.key[ ,.(position = max(position)), by = scaffold] # curious how long the scaffolds are

snv.mat <- as.matrix(snvs[ ,-c(1:4)], rownames = "snv.num")
snv.mat <- t(snv.mat)

snv.dist.bray <- vegdist(x = snv.mat, method = "bray")
snv.nmds.bray <- metaMDS(comm = snv.dist.bray, trymax = 100)
stressplot(snv.nmds.bray)

snv.dist.man <- vegdist(x = snv.mat, method = "manhattan")
snv.nmds.man <- metaMDS(comm = snv.dist.man, trymax = 100)
stressplot(snv.nmds.man)

snv.dist.euc <- vegdist(x = snv.mat, method = "euclidean")
snv.nmds.euc <- metaMDS(comm = snv.dist.euc, trymax = 100)
stressplot(snv.nmds.euc)

dates.key <- unique(snv.instrain[ ,.(sample, date, year, yday, season, invasion)])

# ---- make plots of long-term change ----

# euclidean:
my.nmds <- data.table("sample" = rownames(snv.nmds.euc$points),"x" = snv.nmds.euc$points[ ,1], "y" = snv.nmds.euc$points[ ,2])
# bray-curtis:
my.nmds <- data.table("sample" = rownames(snv.nmds.bray$points),"x" = snv.nmds.bray$points[ ,1], "y" = snv.nmds.bray$points[ ,2])
# manhattan:
my.nmds <- data.table("sample" = rownames(snv.nmds.man$points),"x" = snv.nmds.man$points[ ,1], "y" = snv.nmds.man$points[ ,2])


my.nmds <- merge(x = my.nmds, y = dates.key, by = "sample")
my.nmds[ ,`:=`(season = factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")),
               invasion = factor(invasion, levels = c("none","spiny","zebra")),
               date = parse_date_time(x = date, orders = "ymd"))]
my.nmds <- my.nmds[order(date)]

# ---- for quick re-plotting ----
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-B_nmds_euc.rds")
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-B_nmds_bray.rds")
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-B_nmds_man.rds")
# 
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-A_nmds_euc.rds")
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-A_nmds_bray.rds")
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-A_nmds_man.rds")
# 
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-C_nmds_euc.rds")
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-C_nmds_bray.rds")
# saveRDS(object = my.nmds, file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-C_nmds_man.rds")

my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-B_nmds_euc.rds")
my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-B_nmds_bray.rds")
my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-B_nmds_man.rds")
my.genome <- "acI-B"

my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-A_nmds_euc.rds")
my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-A_nmds_bray.rds")
my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-A_nmds_man.rds")
my.genome <- "acI-A"

my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-C_nmds_euc.rds")
my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-C_nmds_bray.rds")
my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-C_nmds_man.rds")
my.genome <- "acI-C"

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

# ---- end ----

