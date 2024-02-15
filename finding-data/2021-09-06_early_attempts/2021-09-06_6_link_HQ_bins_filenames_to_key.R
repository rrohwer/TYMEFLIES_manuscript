# RRR

stupid <- read.csv(file = "finding-data/2021-09-06_5_HQ_bins_filenames_from_Neha.txt", header = F)
colnames(stupid) <- "Filename.From.Neha"
str(stupid)

key <- readRDS(file = "data/2021-05-24_metaG_stats/all_img_assemblies.rds")
key <- key[ ,c(1,8)]
str(key) 
key <- key[!is.na(key$IMG.Taxon.ID), ]
key <- unique(key)

binID <- sub(pattern = "\\.fna", replacement = "", x = stupid$Filename.From.Neha)

taxonID <- sub(pattern = "_.*\\.fna", replacement = "", x = stupid$Filename.From.Neha)

binnum <- sub(pattern = "^.*_", replacement = "", x = binID)

stupid <- cbind(stupid, "Taxon.ID" = taxonID, "Bin.ID" = binID, "Bin.Number" = binnum)

less.stupid <- merge(x = stupid, y = key, by.x = "Taxon.ID", by.y = "IMG.Taxon.ID", all.x = TRUE, all.y = FALSE)

less.stupid <- less.stupid[ ,c(2,3,1,4,5)]

less.stupid <- cbind(less.stupid, "New.Filename" = paste0(less.stupid$TYMEFLIES.Name,"_",less.stupid$Bin.Number,".fasta"))

# write.csv(x = less.stupid, file = "finding-data/2021-09-06_7_FILENAME_KEY_HQbins_from_Neha.csv", quote = F, row.names = F)
