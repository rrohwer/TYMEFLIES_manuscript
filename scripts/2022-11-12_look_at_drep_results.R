# RRR

folder.drep <- "data/2022-11-10_dRep_files/output/full_run_output/data_tables"
folder.bin.info <- "data/2022-09-23_bin_stats"
folder.output <- "data/2022-11-10_dRep_files/"

# ---- get table of drep results ----

all.genomes <- read.csv(file = file.path(folder.drep, "genomeInformation.csv"))

good.genomes <- read.csv(file = file.path(folder.drep, "Cdb.csv"))
good.genomes <- good.genomes[ ,c("genome", "secondary_cluster")]

all.genomes <- merge(x = all.genomes, y = good.genomes, by = "genome", all = T)

win.genomes <- read.csv(file = file.path(folder.drep, "Wdb.csv"))
win.genomes$winner <- TRUE
win.genomes <- win.genomes[ ,c("genome","score","winner")]

all.genomes <- merge(x = all.genomes, y = win.genomes, by = "genome", all = T)

index <- is.na(all.genomes$winner)
all.genomes$winner[index] <- FALSE

num.in.cluster <- aggregate(rep.int(1, times = nrow(all.genomes)), by = list(all.genomes$secondary_cluster), FUN = sum)
colnames(num.in.cluster) <- c("cluster","num.in.cluster")

all.genomes <- merge(x = all.genomes, y = num.in.cluster, by.x = "secondary_cluster", by.y = "cluster", all = T)

head(all.genomes)

# ---- add gtdb info ----

gtdb <- readRDS(file = file.path(folder.bin.info,"all_bins.rds"))

head(gtdb)

gtdb$Name <- paste0(gtdb$Name,".fna")

gtdb <- gtdb[ ,-(2:3)]

all.genomes <- merge(x = all.genomes, y = gtdb, by.x = "genome", by.y = "Name", all = T)

index <- order(all.genomes$domain, all.genomes$phylum, all.genomes$class, all.genomes$order, all.genomes$family, all.genomes$genus, all.genomes$species, all.genomes$secondary_cluster)
all.genomes <- all.genomes[index, ]


# make a format with just unique ones

winners.only <- all.genomes[all.genomes$winner, ]

# ---- save data ----

saveRDS(object = all.genomes, file = file.path(folder.output,"drep_results_all_genomes.rds"))
saveRDS(object = winners.only, file = file.path(folder.output, "drep_results_winning_genomes.rds"))
write.csv(x = winners.only, file = file.path(folder.output, "drep_results_winning_genomes.csv"), quote = F, row.names = F)
