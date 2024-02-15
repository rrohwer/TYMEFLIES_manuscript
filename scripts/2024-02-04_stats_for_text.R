# RRR

library(data.table)

drep <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv")
drep
median(drep$completeness)
median(drep$contamination)

drep[order == "Nanopelagicales", .N, by = order]
drep[family == "Nanopelagicaceae", .N, by = family]
