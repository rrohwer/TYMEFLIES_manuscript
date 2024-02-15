# RRR

# ---- what's in each file? ----

# 311 genomes were input to dRep
# 111 genomes passed checkM filtering
# 84 genomes were "winners"

my.folder <- "data/2022-11-10_dRep_files/output/test_output/data_tables/"

x <- read.csv(file = file.path(my.folder, "Bdb.csv")) # genome, file location, length
head(x)
nrow(x) # 111

x <- read.csv(file = file.path(my.folder, "Cdb.csv")) # genome, secondary_cluster, threshold, cluster_method ,comparison_algorithm, primary_cluster
head(x)
nrow(x) # 111
x$secondary_cluster # which cluster it's in
# seems like this is a key file, inlcudes each genome and the cluster it belongs to
cbind(x$primary_cluster, x$secondary_cluster) # secondary cluster is the finer-resolution cluster

x <- read.csv(file = file.path(my.folder, "Mdb.csv")) # genome1, genome2, distance, similarity
head(x)
nrow(x) # 12321

x <- read.csv(file = file.path(my.folder, "Ndb.csv")) # reference, querry, ani, alignment_coverage, primary_cluster
head(x)
nrow(x) # 177

x <- read.csv(file = file.path(my.folder, "Sdb.csv")) # genome, score
head(x)
nrow(x) # 111

x <- read.csv(file = file.path(my.folder, "Wdb.csv")) # genome, cluster, score
head(x)
nrow(x) # 84  ONLY WINNERS
x$cluster
# seems promising but is only the winning genomes

x <- read.csv(file = file.path(my.folder, "Widb.csv")) # checkM info, so empty since I provide it separately
head(x)
nrow(x) # 84 ONLY WINNERS
x$cluster

x <- read.csv(file = file.path(my.folder, "genomeInfo.csv")) # genome, completeness, contamination, length, N50
head(x)
nrow(x) # 311 

x <- read.csv(file = file.path(my.folder, "genomeInformation.csv")) # genome, completeness, contamination, length, N50, centrality
head(x)
nrow(x) # 311 
# Completeness, Contamination, and strain heterogeneity are provided by the user or calculated with checkM. 
# N50 is a measure of how big the pieces are that make up the genome. size is the total length of the genome. 
# Centrality is a measure of how similar a genome is to all other genomes in it’s cluster. This metric helps 
# pick genome that are similar to all other genomes, and avoid picking genomes that are relative outliers.


# ---- Try pulling out just the cluster info I want: ----

all.genomes <- read.csv(file = file.path(my.folder, "genomeInformation.csv"))

good.genomes <- read.csv(file = file.path(my.folder, "Cdb.csv"))
good.genomes <- good.genomes[ ,c("genome", "secondary_cluster")]

all.genomes <- merge(x = all.genomes, y = good.genomes, by = "genome", all = T)

win.genomes <- read.csv(file = file.path(my.folder, "Wdb.csv"))
win.genomes$winner <- TRUE
win.genomes <- win.genomes[ ,c("genome","score","winner")]

all.genomes <- merge(x = all.genomes, y = win.genomes, by = "genome", all = T)

head(all.genomes)

# OK I think this will be fine