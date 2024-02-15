x <- read.table(file = "data/2022-11-17_contig_counts/winnernumcontigs.txt", header = F, sep = ":")
y <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")

z <- merge(x = x, y = y, by.x = "V1", by.y = "bin.filename", all.x = T, all.y = F)

finishable <- z[z$completeness > 99, ]
colnames(finishable)[2] <- "num.contigs"
finishable <- finishable[ ,c("bin.full.name","drep.cluster.95ANI","num.in.cluster","completeness","contamination","num.contigs","length","N50","phylum","class","order","family","genus","species")]

write.csv(x = finishable, file = "data/2022-11-17_contig_counts/winners_with_99_completeness.csv", row.names = F, quote = F)

summary(finishable$num.contigs)

finishable <- finishable[finishable$num.contigs < 20, ]

write.csv(x = finishable, file = "data/2022-11-17_contig_counts/winners_with_99_completeness_and_20_contigs.csv", row.names = F, quote = F)

finishable <- finishable[finishable$num.contigs < 5, ]
