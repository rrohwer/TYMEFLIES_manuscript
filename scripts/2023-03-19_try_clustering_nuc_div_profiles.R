# RRR
# Idea is to classify into seasonal patterns, 
# then look for trends in taxonomy or not
# for now focus on seasonal signals
# later can look for long-term signals also 
# clustering could simplify looking through them at least?

# ---- set up ----

library(data.table)
library(ggplot2)
library(lubridate)

genome.info <- readRDS(file = "data/2023-02-12_instrain_results_TIMEOUT/genome_info_combined.rds")
coverm <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")

# ---- try clustering coverm abundance ----

colnames(coverm)[3] <- "Rel.Abund.perc"
unmapped <- coverm[Genome == "unmapped", .(Sample, Rel.Abund.perc)] # same info about percent mapped !

# normalize each bin by its max abundance
abund <- dcast(coverm, Genome ~ Sample, value.var = "Rel.Abund.perc")
abund[1:5,1:5]
row.samples <- abund$Genome
abund <- as.matrix(abund[ ,-1])
row.names(abund) <- row.samples
max.abunds <- rowSums(abund)
abund.perc.max <- abund / max.abunds
rowSums(abund.perc.max)
colSums(abund.perc.max)
abund.perc.max[1:5,1:5]

# try simplest form of kmeans clustering ----
k.clusts <- kmeans(x = abund.perc.max, centers = 100, iter.max = 100, nstart = 1, algorithm = "Lloyd")
k.clusts$cluster[1:5] # which cluster each contig belongs to
head(k.clusts$centers) # matrix of cluster means- the average value for each cluster on each day
k.clusts$size # number of datapoints in each cluster
hist(k.clusts$size)

# parse into a long output table
abund.k.clusts <- data.table("Genome" = row.names(abund.perc.max), abund.perc.max)
abund.k.clusts <- melt(data = abund.k.clusts, id.vars = "Genome")
colnames(abund.k.clusts) <- c("Genome","Sample","Abund.perc.max")
abund.k.clusts$Date <- substr(x = abund.k.clusts$Sample, start = 3, stop = 12) |>
  parse_date_time(orders = "ymd")
clust.key <- data.table("Genome" = names(k.clusts$cluster),"Kcluster" = k.clusts$cluster)
abund.k.clusts <- merge(x = clust.key, y = abund.k.clusts, by = "Genome")

my.clust <- unique(abund.k.clusts$Kcluster)[4]
ggplot(data = abund.k.clusts[Kcluster == my.clust, ], aes(x = Date, y = Abund.perc.max))+
  geom_line(aes(group = Genome))


# try hierarchical clustering ----
abund.dist <- dist(x = abund.perc.max, method = "manhattan")
str(abund.dist) # genomes x genomes distance matrix
h.clusts <- hclust(d = abund.dist, method = "ward.D")
print(h.clusts)
names(h.clusts)
head(h.clusts$labels) # the genome names
head(h.clusts$order) # the order from left to write so branches don't cross
head(h.clusts$height) # the height of the branches
head(h.clusts$merge) # maybe the acutal data about where branches happen?
# to assign genomes to a cluster, need to choose a height cutoff
plot(h.clusts, labels = F)
hist(h.clusts$height)
h.clusts.cut <- cutree(tree = h.clusts, h = 1.5) 
str(h.clusts.cut) # just like the kmean $cluster object
length(unique(h.clusts.cut)) 

# parse into a long output table
abund.h.clusts <- data.table("Genome" = row.names(abund.perc.max), abund.perc.max)
abund.h.clusts <- melt(data = abund.h.clusts, id.vars = "Genome")
colnames(abund.h.clusts) <- c("Genome","Sample","Abund.perc.max")
abund.h.clusts$Date <- substr(x = abund.h.clusts$Sample, start = 3, stop = 12) |>
  parse_date_time(orders = "ymd")
clust.key <- data.table("Genome" = names(h.clusts.cut),"Hcluster" = h.clusts.cut)
abund.h.clusts <- merge(x = clust.key, y = abund.h.clusts, by = "Genome")

my.clust <- unique(abund.h.clusts$Hcluster)[5]
ggplot(data = abund.h.clusts[Hcluster == my.clust, ], aes(x = Date, y = Abund.perc.max))+
  geom_line(aes(group = Genome))

# ---- try clustering nucleotide diversity ----

# normalize each genome by max diversity
genome.info <- data.table(genome.info)
div <- dcast(genome.info, genome ~ sample, value.var = "nucl_diversity")
div[1:5,1:5]
row.samples <- div$genome
div <- as.matrix(div[ ,-1])
row.names(div) <- row.samples
max.divs <- rowSums(div, na.rm = T)
div.perc.max <- div / max.divs
rowSums(div.perc.max, na.rm = T)
colSums(div.perc.max, na.rm = T)
div.perc.max[1:5,1:5]
# presence of NAs doesn't work in hclust, what to do? make them zero? group by season?
div.perc.max[is.na(div.perc.max)] <- 0 # for now, even though not best solution

# try hierarchical clustering ----
div.dist <- dist(x = div.perc.max, method = "manhattan")
str(div.dist) # genomes x genomes distance matrix
h.clusts <- hclust(d = div.dist, method = "ward.D")
# to assign genomes to a cluster, need to choose a height cutoff
plot(h.clusts, labels = F)
abline(h = 2.5, col = "red")
hist(h.clusts$height)
h.clusts.cut <- cutree(tree = h.clusts, h = 2.5) 
str(h.clusts.cut) # just like the kmean $cluster object
length(unique(h.clusts.cut)) 

# parse into a long output table
div.h.clusts <- data.table("Genome" = row.names(div.perc.max), div.perc.max)
div.h.clusts <- melt(data = div.h.clusts, id.vars = "Genome")
colnames(div.h.clusts) <- c("Genome","Sample","Abund.perc.max")
div.h.clusts$Date <- substr(x = div.h.clusts$Sample, start = 3, stop = 12) |>
  parse_date_time(orders = "ymd")
clust.key <- data.table("Genome" = names(h.clusts.cut),"Hcluster" = h.clusts.cut)
div.h.clusts <- merge(x = clust.key, y = div.h.clusts, by = "Genome")

my.clust <- unique(div.h.clusts$Hcluster)[3]
ggplot(data = div.h.clusts[Hcluster == my.clust, ], aes(x = Date, y = Abund.perc.max))+
  geom_line(aes(group = Genome))
