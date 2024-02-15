# RRR
# not really make the dups unique, but choose the 1 unique one per dup and make the unique file

library(readxl)

manual.choices <- read_excel(path = "finding-data/2021-12-23_filtered_reads/dupSPIDS-manual_choices.xlsx")

bigtyme <- read_excel(path = "data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx", sheet = 1, na = "NA")

# ----

bigtyme <- bigtyme[is.na(bigtyme$`Filenames_filter-METAGENOME`), ]
bigtyme$roworder <- 1:nrow(bigtyme)

unique(manual.choices$USE)
index <- which(manual.choices$USE == "YES")
new.unique.file <- manual.choices[index, 1:2]

colnames(bigtyme)
colnames(new.unique.file)
biggertyme <- merge(x = bigtyme, y = new.unique.file, by.x = "ITS_SPID", by.y = "SPID")
biggertyme <- biggertyme[order(biggertyme$roworder), ]

# ----

write.csv(x = biggertyme, file = "data/2021-12-13_IMG_IDs_for_504350/old/2022-01-20_dupSPID_filteredread_filenames.csv", row.names = F, quote = F)

write.table(x = new.unique.file, file = "finding-data/2021-12-23_filtered_reads/dupSPIDs_made_unique.txt", sep = " ", quote = F, row.names = F, col.names = T)
