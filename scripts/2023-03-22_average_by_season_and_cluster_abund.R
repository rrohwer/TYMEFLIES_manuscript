# RRR
# cluster abundance into seasonal and long-term patterns

library(data.table)
library(lubridate)

rel.abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz", sep = "\t")
genome.info <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
bin.info <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds") 
bin.info <- as.data.table(bin.info)

output.season.avs <- "data/2023-03-21_clustering_analysis_results/abund_seasonal_table.tsv.gz"
output.season.all.data <- "data/2023-03-21_clustering_analysis_results/abund_seasonal_table_big.tsv.gz"
output.annual.avs <- "data/2023-03-21_clustering_analysis_results/abund_annual_table.tsv.gz"
output.annual.all.data <- "data/2023-03-21_clustering_analysis_results/abund_annual_table_big.tsv.gz"

plot.folder <-"figures/2023-03-21_visualize_clusters/"

# ---- combine breadth info with abund info ----

# subset just relative abundance from coverM output
colnames(rel.abund)[3] <- "abund.perc"
rel.abund <- rel.abund[Genome != "unmapped", .(sample = Sample, genome = Genome, abund.perc)]

# subset inStrain output to genomes with adequate breadth (>= 50% of expected bases covered)
colnames(genome.info)
genome.info <- genome.info[ , .(genome, sample, date, year, yday, season, breadth, breadth_expected, breadth_minCov)]
hist(genome.info$breadth / genome.info$breadth_expected)
genome.info <- genome.info[ (breadth / breadth_expected) >= .5, ] # **** is this too lenient?? ****
genome.info <- genome.info[ ,-c("breadth", "breadth_expected", "breadth_minCov")]

# combine and subset abund to only genomes with enough breadth
rel.abund <- merge(x = genome.info, y = rel.abund, by = c("genome","sample"), all.x = TRUE, all.y = FALSE)
rel.abund

# ---- average abundance by season ----

# average for each season in each year
by.season <- rel.abund[order(genome, year, season), .(abund.perc = mean(abund.perc, na.rm = T)), by = .(genome, year, season)]

# average seasons across all years
by.season.overlay <- by.season[ , .(abund.perc = mean(abund.perc, na.rm = T)), by = .(genome, season)]

# average years across all seasons, and normalized by max value
by.year.averages <- by.season[ , .(abund.perc = mean(abund.perc, na.rm = T)), by = .(genome, year)]


# ---- convert to matrices ----

# make matrix
seasonal <- dcast(data = by.season.overlay, formula = genome ~ season, value.var = "abund.perc")
seasonal <- as.matrix(seasonal, rownames = T)

annual <- dcast(data = by.year.averages, formula = genome ~ year, value.var = "abund.perc")
annual <- as.matrix(annual, rownames = T)

# normalize by genome max values
my.maxs <- apply(seasonal, 1, max, na.rm = T)
seasonal <- seasonal / my.maxs * 100

my.maxs <- apply(annual, 1, max, na.rm = T)
annual <- annual / my.maxs * 100

# include only genomes present in all seasons / all years 
index <- apply(is.na(seasonal), 1, any)
seasonal <- seasonal[!index, ]
index <- apply(seasonal == 0, 1, any) # none? all zeros must be NA
seasonal <- seasonal[!index, ]

index <- which(colnames(annual) == "2019") # exclude 2019 partial year
annual <- annual[, -index]
index <- apply(is.na(annual), 1, any)
annual <- annual[!index, ]
index <- apply(annual == 0, 1, any) # none? all zeros must be NA
annual <- annual[!index, ]


# ---- perform hierarchical clustering ----

# get distance matrix
seasonal.dist <- dist(x = seasonal, method = "manhattan")

annual.dist <- dist(x = annual, method = "manhattan")

# get clustering results
seasonal.clust <- hclust(d = seasonal.dist, method = "ward.D")

annual.clust <- hclust(d = annual.dist, method = "ward.D")

# choose cluster height cutoffs 
# 1500
# cutoff = 1500 # gives 6, but don't correspond to seasons exactly
cutoff = 725 # looks pretty good 
plot(seasonal.clust, labels = F)
abline(h = cutoff, col = "red")
text(x = 0, y = cutoff, labels = length(unique(cutree(seasonal.clust, h = cutoff))))
seasonal.clust.cut <- cutree(tree = seasonal.clust, h = cutoff)
seasonal.clust.cut <- data.table(genome = names(seasonal.clust.cut), cluster = seasonal.clust.cut)

# cutoff = 600 # looks pretty good, gives 10 clusters
cutoff = 1000 # better defined, can see the big-change ones better
plot(annual.clust, labels = F)
abline(h = cutoff, col = "red")
text(x = 0, y = cutoff, labels = length(unique(cutree(annual.clust, h = cutoff))))
annual.clust.cut <- cutree(tree = annual.clust, h = cutoff)
annual.clust.cut <- data.table(genome = names(annual.clust.cut), cluster = annual.clust.cut)

# expand tables of cluster membership
seasonal.table <- data.table(genome = row.names(seasonal), seasonal)
seasonal.table <- merge(x = seasonal.clust.cut, y = seasonal.table, by = "genome", all = T)
seasonal.table <- merge(x = seasonal.table, y = bin.info, by.x = "genome", by.y = "bin.full.name", all.x = TRUE, all.y = FALSE)
seasonal.table <- melt(data = seasonal.table, measure.vars = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), 
                       variable.name = "season", value.name = "abund.perc.max")
seasonal.table <- merge(x = seasonal.table, y = by.season.overlay, by = c("genome","season"))

annual.table <- data.table(genome = row.names(annual), annual)
annual.table <- merge(x = annual.clust.cut, y = annual.table, by = "genome", all = T)
annual.table <- merge(x = annual.table, y = bin.info, by.x = "genome", by.y = "bin.full.name", all.x = TRUE, all.y = FALSE)
annual.table <- melt(data = annual.table, measure.vars = c("2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018"), 
                     variable.name = "year", value.name = "abund.perc.max")
annual.table[ , year := as.numeric(as.character(year))]
annual.table <- merge(x = annual.table, y = by.year.averages, by = c("genome","year"))

# make tables with raw data values and cluster membership
seasonal.big.table <- rel.abund[order(genome, year, season), .(genome, abund.perc, date, year, yday, season)]
seasonal.big.table <- merge(x = seasonal.clust.cut, y = seasonal.big.table, by = "genome", all.x = TRUE, all.y = FALSE)
seasonal.big.table <- merge(x = seasonal.big.table, y = bin.info, by.x = "genome", by.y = "bin.full.name", all.x = TRUE, all.y = FALSE)

annual.big.table <- rel.abund[order(genome, year, season), .(genome, abund.perc, date, year, yday, season)]
annual.big.table <- merge(x = annual.clust.cut, y = annual.big.table, by = "genome", all.x = TRUE, all.y = FALSE)
annual.big.table <- merge(x = annual.big.table, y = bin.info, by.x = "genome", by.y = "bin.full.name", all.x = TRUE, all.y = FALSE)

# ---- save results files for plotting ----

fwrite(x = seasonal.table, file = output.season.avs, sep = "\t", compress = "gzip")
fwrite(x = seasonal.big.table, file = output.season.all.data, sep = "\t", compress = "gzip")

fwrite(x = annual.table, file = output.annual.avs, sep = "\t", compress = "gzip")
fwrite(x = annual.big.table, file = output.annual.all.data, sep = "\t", compress = "gzip")


# ---- OK but some quick looks now ----


library(ggplot2)

seasonal.table[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
seasonal.table[ ,season.num := as.numeric(season)]
col.key <- seasonal.table[ , .N, by = phylum]
col.key <- col.key[order(-N)]
col.key$color <- c("darkred","blue2","orange","purple3","green3","grey","magenta3","cyan2","grey","grey","grey","grey","grey")
col.key[is.na(phylum), phylum := "NA"]
seasonal.table[is.na(phylum), phylum := "NA"]
seasonal.table[ ,phylum := factor(phylum, levels = paste0(col.key[ ,phylum]))]

p1 <- ggplot(data = seasonal.table, aes(x = season.num, y = abund.perc.max, color = phylum))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = col.key$color)+
  geom_vline(xintercept = 1:6, alpha = .25)+
  geom_line(aes(group = genome), alpha = .5)+
  facet_wrap(~cluster)+
  scale_x_continuous(name = element_blank(), breaks = 1:6, labels = levels(seasonal.table$season))+
  guides(color = guide_legend(override.aes = list(linewidth = 5, alpha = 1)))

p1a <- ggplot(data = seasonal.table, aes(x = season.num, y = abund.perc, color = phylum))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = col.key$color)+
  geom_vline(xintercept = 1:6, alpha = .25)+
  geom_line(aes(group = genome), alpha = .5)+
  facet_wrap(~cluster)+
  scale_x_continuous(name = element_blank(), breaks = 1:6, labels = levels(seasonal.table$season))+
  guides(color = guide_legend(override.aes = list(linewidth = 5, alpha = 1)))


annual.key <- annual.table[ , .N, by = phylum]
annual.key <- merge(x = col.key, y = annual.key, by = "phylum", all.x = F, all.y = T)
annual.key[ ,N.x := NULL]
annual.key <- annual.key[order(-N.y)]
annual.key[is.na(phylum), phylum := "NA"]
annual.table[is.na(phylum), phylum := "NA"]
annual.table[ ,phylum := factor(phylum, levels = annual.key[ ,phylum])]

p2 <- ggplot(data = annual.table, aes(x = year, y = abund.perc.max, color = phylum))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = annual.key$color)+
  # geom_vline(xintercept = 2000:2018, alpha = .25)+
  geom_line(aes(group = genome), alpha = .5)+
  facet_wrap(~cluster)+
  scale_x_continuous(name = element_blank(), breaks = 2000:2018)+
  guides(color = guide_legend(override.aes = list(linewidth = 5, alpha = 1)))

p2a <- ggplot(data = annual.table, aes(x = year, y = abund.perc, color = phylum))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = annual.key$color)+
  # geom_vline(xintercept = 2000:2018, alpha = .25)+
  geom_line(aes(group = genome), alpha = .5)+
  facet_wrap(~cluster)+
  scale_x_continuous(name = element_blank(), breaks = 2000:2018)+
  guides(color = guide_legend(override.aes = list(linewidth = 5, alpha = 1)))

png(filename = file.path(plot.folder,"abund_seasonal_clusts-normalized_div.png"), width = 12, height = 9, units = "in", res = 300)
p1
dev.off()

png(filename = file.path(plot.folder,"abund_seasonal_clusts-true_div.png"), width = 12, height = 9, units = "in", res = 300)
p1a
dev.off()

png(filename = file.path(plot.folder,"abund_annual_clusts-normalized_div.png"), width = 12, height = 9, units = "in", res = 300)
p2
dev.off()

png(filename = file.path(plot.folder,"abund_annual_clusts-true_div.png"), width = 12, height = 9, units = "in", res = 300)
p2a
dev.off()
