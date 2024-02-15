# RRR
# for unique IDS, here's where the filenames got added to bigtyme
# but for dup spids, already added the filenames in the script that made them unique

# Still working on getting & organizing the filtered read metagenome files

library(readxl)

# file.uniqueSPIDs <- "finding-data/2021-12-23_filtered_reads/uniqueSPIDS.txt"
file.uniqueSPIDs <- "finding-data/2021-12-23_filtered_reads/dupSPIDs_made_unique.txt"
file.bigtyme <- "data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx"

# file.biggertyme <- "data/2021-12-13_IMG_IDs_for_504350/old/2022-01-05_bigtyme_with_reads_filenames-uniqueSPIDs.csv"
file.newscript <- "finding-data/2021-12-23_filtered_reads/2022-01-20_6_rename_files-dupSPIDs.sh"

# ----

bigtyme <- read_excel(path = file.bigtyme, sheet = 1, na = "NA")
bigtyme <- as.data.frame(bigtyme)
bigtyme$row.order <- c(1:nrow(bigtyme))
head(bigtyme)

spids <- read.delim(file = file.uniqueSPIDs, header = T, sep = " ")
head(spids)


bigtyme <- merge(x = spids, y = bigtyme, by.x = "SPID", by.y = "ITS_SPID", all = T)
head(bigtyme)
bigtyme <- bigtyme[order(bigtyme$row.order), ]

bigtyme$new.names <- paste(bigtyme$TYMEFLIES.name,
                            bigtyme$taxon_oid,
                            "filtered-reads.fastq.gz",
                           sep = "_")

bigtyme.trim <- bigtyme[!is.na(bigtyme$FILENAME), ]
bigtyme.trim$FILENAME <- sub(pattern = "^.*/", replacement = "", x = bigtyme.trim$FILENAME)

new.script <- paste("mv", bigtyme.trim$FILENAME, bigtyme.trim$new.names)

# ---- export ----

# write.csv(x = bigtyme, file = file.biggertyme, quote = F, row.names = F)
write.table(x = new.script, file = file.newscript, quote = F, row.names = F, col.names = F)
