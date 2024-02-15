# RRR
# exported summary of # genes per COG category that are under selection on diff dates
# have one summary file per genome
# combine all the genomes into the same big file

# abund from coverM is the only data object that has everything. 
# InStrain leaves samples out if they're below a threshold.
# so need to make table to that backbone to not have blanks

library(data.table)

# ---- import ----

# get tax
tax <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
tax <- as.data.table(tax)
tax <- tax[winner == TRUE, .(bin.full.name, domain, phylum, class, order, family, genus, species)]

# get COGs
folder.path <- "data/2023-05-31_COG_summary_data/Selection_Summaries/"
my.files <- list.files(folder.path)

all.genomes.list <- list()
for (f in my.files){
  one.g <- fread(file = file.path(folder.path,f))
  if (nrow(one.g > 0)){ # some had no selection so are empty tables
    one.g$genome <- sub(pattern = "_selected_COGs\\.tsv\\.gz", replacement = "", x = f)
    all.genomes.list <- c(all.genomes.list, list(one.g))
  }
}

all.genomes.table <- rbindlist(l = all.genomes.list)
all.genomes.table <- all.genomes.table[ ,-c("date", "season", "year")]

# get abunds
abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz", sep = "\t")
colnames(abund)[3] <- "rel.abund.perc"
abund <- abund[Genome != "unmapped", .(Genome, Sample,rel.abund.perc)]
 
# get breadth and nucl div
breadth <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
breadth <- breadth[ ,.(genome, breadth, breadth_expected, nucl_diversity, sample, date, year, yday, season)]

# ---- combine ----

combo <- merge(x = abund, y = breadth, by.x = c("Sample", "Genome"), by.y = c("sample","genome"), all.x = T, all.y = T)
combo <- merge(x = combo, y = tax, by.x = c("Genome"), by.y = c("bin.full.name"), all.x = T, all.y = T)
combo <- merge(x = combo, y = all.genomes.table, by.x = c("Sample", "Genome"), by.y = c("sample","genome"))
