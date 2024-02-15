

library(readxl)

file.chosenfiles <- "finding-data/2022-03-08_assembled_contigs/dupSPIDS_chosen.xlsx"

created.file <- "finding-data/2022-03-08_assembled_contigs/dupSPIDS_chosen.txt"


dups <- read_excel(path = file.chosenfiles, n_max = 12)

index <- which(dups$USE == "YES")

dups <- dups[index, ]

dups <- dups[ ,1:2]

write.table(x = dups, file = created.file, quote = F, sep = " ", row.names = F, col.names = T)
