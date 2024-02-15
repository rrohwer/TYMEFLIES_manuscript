# RRR
# similar to 2022-11-12_look_at_drep_results.R

folder.drep <- "data/2023-02-24_dRep_with_range_ANIs/output_raw/dRep_ANIrange_output"
file.gtdb <- "data/2022-09-23_bin_stats/all_bins.rds"
folder.output <- "data/2023-02-24_dRep_with_range_ANIs/output_processed"

# ---- get table of drep results ----

ani.range <- seq(from = .90, to = .99, by = .01)

for (ani in ani.range){
  
  # ---- info from dRep
  
  genome.file <- paste0("genomeInformation_",ani,".csv")
  wdb.file <- paste0("Wdb_",ani,".csv")
  cdb.file <- paste0("Cdb_",ani,".csv")
  
  all.genomes <- read.csv(file = file.path(folder.drep, genome.file))
  
  good.genomes <- read.csv(file = file.path(folder.drep, cdb.file))
  good.genomes <- good.genomes[ ,c("genome", "secondary_cluster")]
  
  all.genomes <- merge(x = all.genomes, y = good.genomes, by = "genome", all = T)
  
  win.genomes <- read.csv(file = file.path(folder.drep, wdb.file))
  win.genomes$winner <- TRUE
  win.genomes <- win.genomes[ ,c("genome","score","winner")]
  
  all.genomes <- merge(x = all.genomes, y = win.genomes, by = "genome", all = T)
  
  index <- is.na(all.genomes$winner)
  all.genomes$winner[index] <- FALSE
  
  num.in.cluster <- aggregate(rep.int(1, times = nrow(all.genomes)), by = list(all.genomes$secondary_cluster), FUN = sum)
  colnames(num.in.cluster) <- c("cluster","num.in.cluster")
  
  all.genomes <- merge(x = all.genomes, y = num.in.cluster, by.x = "secondary_cluster", by.y = "cluster", all = T)
  
  # ---- info from gtdb-tk
  
  gtdb <- readRDS(file = file.gtdb)
  
  gtdb$Name <- paste0(gtdb$Name,".fna")
  
  gtdb <- gtdb[ ,-(2:3)]
  
  all.genomes <- merge(x = all.genomes, y = gtdb, by.x = "genome", by.y = "Name", all = T)
  
  index <- order(all.genomes$domain, all.genomes$phylum, all.genomes$class, all.genomes$order, all.genomes$family, all.genomes$genus, all.genomes$species, all.genomes$secondary_cluster)
  all.genomes <- all.genomes[index, ]
  
  winners.only <- all.genomes[all.genomes$winner, ]
  
  # ---- add full key info
  
  colnames(all.genomes)[1] <- "bin.filename"
  colnames(all.genomes)[2] <- paste0("drep.cluster.",ani,"ANI")
  
  all.genomes$bin.full.name <- sub(pattern = "\\.fna$", replacement = "", x = all.genomes$bin.filename)
  
  bin.names <- strsplit(x = all.genomes$bin.full.name, split = "_")
  x <- t(matrix(unlist(bin.names), nrow=4))
  all.genomes$sample.name <- x[,1]
  all.genomes$taxon.id <- x[,2]
  all.genomes$tymeflies.name <- paste0(x[,1],"_",x[,2])
  all.genomes$binning.group <- x[,3]
  all.genomes$bin.id <- x[,4]
  
  cbind(1:ncol(all.genomes), colnames(all.genomes))
  index <- c(21,19,20,22,23,18,1,3,4,2,5:17)
  cbind(1:ncol(all.genomes), colnames(all.genomes)[index])
  all.genomes <- all.genomes[ ,index]
  head(all.genomes)
  
  # ---- save data
  
  saveRDS(object = all.genomes, file = file.path(folder.output, paste0("drep_results_all_genomes_",ani,".rds")))
  write.csv(x = winners.only, file = file.path(folder.output, paste0("drep_results_winning_genomes_",ani,".csv")), quote = F, row.names = F)
}



