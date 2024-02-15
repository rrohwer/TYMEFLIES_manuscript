# RRR

bins <- readRDS("data/2022-11-10_dRep_files/drep_results_all_genomes.rds")
head(bins)

colnames(bins)[1] <- "bin.filename"
colnames(bins)[2] <- "drep.cluster.95ANI"

bins$bin.full.name <- sub(pattern = "\\.fna$", replacement = "", x = bins$bin.filename)
bins$sample.name <- NA
bins$taxon.id <- NA
bins$tymeflies.name <- NA
bins$binning.group <- NA
bins$bin.id <- NA

bin.names <- strsplit(x = bins$bin.full.name, split = "_")

# <<< SPEED IMPROVEMENTS BY "FAT WILD TUNA" >>>
# --::^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^::--
x <- t(matrix(unlist(bin.names), nrow=4))
head(x)
bins$sample.name <- x[,1]
bins$taxon.id <- x[,2]
bins$tymeflies.name <- paste0(x[,1],"_",x[,2])
bins$binning.group <- x[,3]
bins$bin.id <- x[,4]
# ~~~~~~~~~~~~~~( goodbye Snail .. )~~~~~~~~~~~~~~

#for (e in 1:length(bin.names)){
#  bins$sample.name[e] <- bin.names[[e]][[1]]
#  bins$taxon.id[e] <- bin.names[[e]][2]
#  bins$tymeflies.name[e] <- paste0(bin.names[[e]][[1]],"_",bin.names[[e]][2])
#  bins$binning.group[e] <- bin.names[[e]][3]
#  bins$bin.id[e] <- bin.names[[e]][4]
#}

head(bins)
cbind(1:ncol(bins), colnames(bins))
cbind(1:ncol(bins), colnames(bins)[c(21,19,20,22,23,18,1,3,4,2,5:17)])
bins <- bins[ ,c(21,19,20,22,23,18,1,3,4,2,5:17)]
head(bins)

saveRDS(object = bins, file = "data/2022-11-15_bin_stats/all_bins.rds")

