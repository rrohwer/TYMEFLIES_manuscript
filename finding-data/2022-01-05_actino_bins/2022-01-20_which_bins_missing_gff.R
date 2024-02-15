
library(magrittr)

gff <- read.delim(file = "data/2022-01-05_IMG_IDs_for_Actino_bins/all_gff_files.txt", header = F)
fasta <- read.delim(file = "data/2022-01-05_IMG_IDs_for_Actino_bins/all_fasta_files.txt", header = F)
with16S <- read.delim(file = "data/2022-01-05_IMG_IDs_for_Actino_bins/all_16S_files.txt", header = F)

gff <- gff$V1 %>%
  sub(pattern = "\\.gff", replacement = "", x = .)

fasta <- fasta$V1 %>%
  sub(pattern = "\\.fna", replacement = "", x = .)

with16S <- with16S$V1 %>%
  sub(pattern = "_16S\\.fasta", replacement = "", x = .)

bins.missing.gff <- setdiff(x = fasta, y = gff)

bins.with.gff <- intersect(x = fasta, y = gff)

bins.missing.16S <- setdiff(x = bins.with.gff, y = with16S)

write.table(x = bins.missing.16S, file = "data/2022-01-05_IMG_IDs_for_Actino_bins/no_16S_found.txt", quote = F, row.names = F, col.names = F)
write.table(x = bins.missing.gff, file = "data/2022-01-05_IMG_IDs_for_Actino_bins/no_gff_found.txt", quote = F, row.names = F, col.names = F)

taxon.no.16S <- bins.missing.16S %>%
  sub(pattern = "_.*$", replacement = "",x = .) %>%
  unique(.)

write.table(x = taxon.no.16S, file = "data/2022-01-05_IMG_IDs_for_Actino_bins/no_gff_found-taxonOID.txt", quote = F, row.names = F, col.names = F)
