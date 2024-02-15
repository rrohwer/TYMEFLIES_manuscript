# RRR
# As far as I can tell, there is no way to search for taxon OID in jamo
# so do SPID with some massive grepping

library(readxl)

bigtyme <- read_excel(path = "data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx", sheet = 1, na = "NA")
binids <- read_excel(path = "data/2022-01-05_IMG_IDs_for_Actino_bins/IMG_bins_Actinobacteria.xlsx")

bigtyme <- as.data.frame(bigtyme)
head(bigtyme)
binids <- as.data.frame(binids)
head(binids)

bigids <- merge(x = binids, y = bigtyme, by.x = "IMG Genome ID", by.y = "taxon_oid", all.x = T)
head(bigids)

index.keep <- which(bigids$`16s rRNA` > 0)
bigids <- bigids[index.keep, ] # only 600 actino bins that have 16S

colnames(bigids)

new.script <- paste("jamo info all spid", bigids$ITS_SPID, "| grep", bigids$`Bin ID`, ">> 1_jamo_info.txt")
new.script[1] <- sub(pattern = ">>", replacement = ">", x = new.script[1])

write.csv(x = bigids, file = "data/2022-01-05_IMG_IDs_for_Actino_bins/IMG_bins_Actinobacteria_with16S.csv")
write.table(x = new.script, file = "finding-data/2022-01-05_actino_bins/2022-01-05_1_jamo_info_bins.sh", quote = F, row.names = F, col.names = F)
