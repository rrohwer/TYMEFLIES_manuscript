# RRR

# first get target reps, plus some extras in case checkM comes back bad
# so get 1/family for the "other" orders
# and 1/species for the "nanopelagicales" order
# will run checkM2 and phylosift on all those, and then trim down to best reps

# get the genbank summary file:
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
# https://ftp.ncbi.nlm.nih.gov/genomes/genbank/

library(data.table)
library(tidyr)

set.seed(85684)

gtdb <- fread(file = "data/2023-04-16_GTDB_accessions/1a-Actinomycetia-gtdb-search.tsv", sep = "\t")

gtdb <- separate(data = gtdb, col = "ncbi_taxonomy", into = c("ncbi.domain","ncbi.phylum","ncbi.class","ncbi.order","ncbi.family","ncbi.genus","ncbi.species"), sep = "; ", remove = T, extra = "warn", fill = "right")
gtdb <- separate(data = gtdb, col = "gtdb_taxonomy", into = c("gtdb.domain","gtdb.phylum","gtdb.class","gtdb.order","gtdb.family","gtdb.genus","gtdb.species"), sep = "; ", remove = T, extra = "warn", fill = "right")

gtdb <- as.data.table(gtdb)

genbank <- fread(file = "data/2023-04-16_GTDB_accessions/1a-assembly_summary_genbank.txt.gz", sep = "\t", na.strings = "na", skip = 1, quote = F)

# remove options that are missing in the genbank file
gtdb <- merge(x = gtdb, y = genbank, by.x = "accession", by.y = "# assembly_accession", all = F) # gtdb goes from 27901 refs to 6451

# for each "nano" order, choose 3 per species if exist (incl the species rep)
nano.rep <- gtdb[gtdb.order == "o__Nanopelagicales" & gtdb_species_representative == TRUE]
nano.extra <- gtdb[gtdb.order == "o__Nanopelagicales"  & gtdb_species_representative == FALSE, .N, by = gtdb.species]
nano.extra <- merge(x = nano.extra, y = gtdb, by = "gtdb.species", all.x = T, all.y = F)
nano.extra <- nano.extra[gtdb_species_representative == FALSE]
nano.few <- nano.extra[N < 2, ]
nano.lots <- nano.extra[N >= 2, .(accession = sample(accession, size = 2, replace = F)), by = gtdb.species]
nano <- rbind(nano.rep[ ,.(accession)], nano.few[ ,.(accession)], nano.lots[ ,.(accession)])

# for each "other" order, narrow down to 1 per family (incl only species reps) (later will chose 2-3 family reps per order)
other <- gtdb[gtdb_species_representative == TRUE & gtdb.order != "o__Nanopelagicales", .(accession = sample(accession, size = 1, replace = F)), by = gtdb.family]
other <- merge(x = gtdb, y = other, by = c("accession","gtdb.family"), all.x = F, all.y = T)
other <- other[ ,.(accession)]

# all accessions to retrieve:
get.these <- rbind(other, nano)

get.these.key <- merge(x = get.these, y = gtdb, by = "accession", all.x = T, all.y = F)

# --- to do ----

fwrite(file = "data/2023-04-16_GTDB_accessions/1b-GTDB_download_Actinomycetia_Key.tsv", x = get.these.key, sep = "\t")

# get these
# run checkM2 on them
# narrow down which to actually use 

# ---- sanity check ----

# if do 3 per other, that's

num.other <- gtdb[gtdb.order != "f__Nanopelagicaceae" & gtdb_species_representative == TRUE, .N, by = gtdb.order]
num.other[N >3, N := 3]
num.other <- sum(num.other[ ,N]) # 121 (107 when only include ones with accessions)

num.nano <- gtdb[gtdb.family == "f__Nanopelagicaceae" & gtdb_species_representative == TRUE, .N, by = gtdb.species]
num.nano <- sum(num.nano[ ,N]) # 216 if do all nano species (205 when only include oens with accessions)

tot.backbone <- num.other + num.nano # 337 (312 when only include ones with accessions)

# and how many MAGs added in?

tax <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
tax <- data.table(tax)
tax <- tax[winner == TRUE & class == "c__Actinomycetia"]
num.tax <- nrow(tax) # 256
