# RRR
# 9/25/22
# combine checkm and gtdb output into one bin stats file

checkm <- readRDS("data/2022-09-17_checkm2_output/combined_checkm2_quality_files.rds")
gtdb <- readRDS(file = "data/2022-09-22_gtdb-tk_output/bin_taxonomy.rds")

output.folder <- "data/2022-09-23_bin_stats/"

head(checkm)
head(gtdb)

checkm <- checkm[ ,1:3]
gtdb <- as.data.frame(gtdb) # some bins didn't have enough markers to make an assignment

colnames(checkm)
colnames(gtdb)
bins <- merge(x = checkm, y = gtdb, by.x = "Name", by.y = "binID", all = TRUE)

saveRDS(object = bins, file = file.path(output.folder,"all_bins.rds"))
write.csv(x = bins, file = file.path(output.folder,"all_bins.csv"), row.names = F, quote = F)
