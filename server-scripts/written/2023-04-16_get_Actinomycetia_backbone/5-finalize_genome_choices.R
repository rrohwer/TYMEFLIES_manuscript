# RRR

library(data.table)
library(readxl)

checkm <- fread(file = "data/2023-04-16_GTDB_Actinomycetia_backbone/4-checkM2_quality_report.tsv")

gtdb <- fread(file = "data/2023-04-16_GTDB_Actinomycetia_backbone/1b-GTDB_download_Actinomycetia_Key.tsv")

# ---- check out these files ----

unique(gtdb$gtdb.class)

# oh shit. I totally did all the orders not just orders in class Actinomycetia
# that is weird though, why would I have to filter to Actinomycetia when I only downloaded Actinomycetia???
# so the download feature obviously didn't work? but then why weren't there MORE?

gtdb <- gtdb[gtdb.class == "c__Actinomycetia"] # OK, so obviously there were just a few extras?

unique(gtdb$gtdb.class)
unique(gtdb$gtdb.order)

# ---- merge checkm results into gtdb key ----

# all accession numbers end in .1
length(gtdb$accession) # 586
length(grep(pattern = "\\.1", x = gtdb$accession)) # 586

checkm$accession <- sub(pattern = "\\.1_.*$", replacement = "", x = checkm$Name)
checkm$accession <- paste0(checkm$accession,".1")

gtdb <- merge(x = checkm, y = gtdb, by = "accession", all.x = F, all.y = T)

# ---- manually choose which reps to use ----

nano <- gtdb[gtdb.order == "o__Nanopelagicales", .(accession, gtdb.species, gtdb_species_representative, Completeness, Contamination)]

fwrite(x = nano, file = "data/2023-04-16_GTDB_Actinomycetia_backbone/5-manually_choose_nano_reps.csv", sep = ",")
# chose manually, 1 per species. deferred to species rep if it was close, otherwise did the better checkM2

other <- gtdb[gtdb.order != "o__Nanopelagicales", .(accession, gtdb.order, gtdb.family, gtdb_species_representative, Completeness, Contamination)]

fwrite(x = other, file = "data/2023-04-16_GTDB_Actinomycetia_backbone/5-manually_choose_order_reps.csv", sep = ",")
# honestly there are only 44 of these, all different familiees. I think I should just do them all, then at least there is some more resolution

# ---- get final list of bin names, to pull into a backbone reference folder ----

my.nano <- read_excel("data/2023-04-16_GTDB_Actinomycetia_backbone/5-manually_choose_nano_reps.xlsx")
my.nano <- as.data.table(my.nano)
my.nano <- my.nano[!is.na(choose)]

keep <- data.table("accession" = c(my.nano$accession, other$accession))

keep.key <- merge(x = keep, y = gtdb, by = "accession", all.x = T, all.y = F)

# ---- make script to copy out those chosen bins ----

script <- c("#!/bin/bash",
            paste0("cp GTDB_refs/",keep.key$Name,".gz Actinomycetia_backbone/"))

# ---- save ----

fwrite(x = keep.key, file = "data/2023-04-16_GTDB_Actinomycetia_backbone/5b_Actinomycetia_Backbone_Key.tsv", sep = "\t")

write.table(x = script, file = "server-scripts/generated/2023-04-16_get_Actinomycetia_backbone/5-copy_chosen_refs.sh", quote = F, row.names = F, col.names = F)
