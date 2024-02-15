# RRR

binIDs <- read.delim(file = "finding-data/2022-01-05_actino_bins/binIDs.txt")
binIDs <- binIDs[ ,1]

new.script <- paste("Rscript 2022-01-06_4_pull_16S_into_new_fasta.R", binIDs, "bin_gffs bin_fastas bin_16S_fastas")

write.table(x = new.script, file = "finding-data/2022-01-05_actino_bins/2022-01-06_4_pull_16S_from_bins.sh", row.names = F, col.names = F, quote = F)