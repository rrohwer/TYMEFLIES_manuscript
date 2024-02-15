# RRR

library(readxl)

# In excel, sort bigtyme by the SPID
bigtyme <- read_excel("data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx")

key.filenames.uniq <- readRDS("finding-data/2022-03-08_assembled_contigs/filename_key_unique.rds")
key.filenames.dup <- readRDS("finding-data/2022-03-08_assembled_contigs/filename_key_dupSPIDs.rds")

key <- rbind(key.filenames.uniq, key.filenames.dup)
key <- key[ ,-3]
head(key)
key$FILENAME

bigtyme <- merge(x = bigtyme, y = key, by = "ITS_SPID", all = TRUE)

colnames(bigtyme)
colnames(bigtyme)[39] <- "Filenames_assembly.contigs"

bigtyme$Filenames_assembly.contigs

index <- order(bigtyme$ITS_SPID)

bigtyme <- bigtyme[index, ]

write.table(x = bigtyme, file = "data/2021-12-13_IMG_IDs_for_504350/old/2022-03-08_bigtyme_plus_assembly_filenames.tsv",
          quote = FALSE, row.names = FALSE, sep = "\t")
