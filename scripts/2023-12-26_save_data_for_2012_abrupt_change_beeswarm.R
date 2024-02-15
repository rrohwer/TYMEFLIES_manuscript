# RRR


library(data.table)
library(lubridate)
library(beeswarm)

# snv.stats <- fread(file = "data/2023-12-04_long-term_change_classifier/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")
# tax <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
# genome.info <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
# abunds <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
# sample.key <- fread("data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv", colClasses = c("date" = "character"))

snv.stats <- fread(file = "data/2023-12-04_long-term_change_classifier/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")
tax.info <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv")
sample.key <- fread(file = "data/2023-12-08_combined_genome_info/sample_key.tsv")

# ---- combine tables ----

sample.key[ ,date := parse_date_time(date, orders = "ymd")]
sample.key <- sample.key[lesser.duplicate == FALSE]

snv.stats <- merge(x = snv.stats, y = tax.info, by = "genome")

snv.stats[ ,`:=`(max.disturbance.start.date = parse_date_time(max.disturbance.start.date, orders = "ymd"),
                 breakpoint.date = parse_date_time(breakpoint.date, orders = "ymd"))]

summary(snv.stats$num.dates)
summary(snv.stats$num.dates.cov10x)
summary(snv.stats$mean.abund)
nrow(snv.stats)

# ---- pull out the year of abrupt changes ----

snv.changes <- snv.stats[Classified.LT.Change == "step" | Classified.LT.Change == "disturbance"]
snv.changes[Classified.LT.Change == "disturbance", change.date := decimal_date(max.disturbance.start.date)]
snv.changes[Classified.LT.Change == "step", change.date := decimal_date(breakpoint.date)]

# ---- get plotting data ----

snvs <- snv.changes[ ,.(genome,domain,phylum,class,order,family,genus,species,Classified.LT.Change,change.date,mean.abund)]

snvs[Classified.LT.Change == "disturbance", color := "orange2"]
snvs[Classified.LT.Change == "step", color := "magenta3"]

snvs[ ,cex := 1 + mean.abund * 2]

# group by phylum but with acI genera split out
snvs[phylum != "Actinobacteriota", "group.ID" := phylum]
snvs[phylum == "Actinobacteriota" & family != "Nanopelagicaceae", "group.ID" := "Other"]
snvs[family == "Nanopelagicaceae", "group.ID" := genus]

bee.list <- list()
for(g in unique(snvs$group.ID)){
  bee <- swarmx(x = 0, y = snvs[group.ID == g, change.date], priority = "random", cex = 2)
  bee <- data.table("genome" = snvs[group.ID == g, genome], "x.shift" = bee$x)
  bee.list <- c(bee.list, list(bee))
}
bee.table <- rbindlist(bee.list)

snvs <- merge(y = bee.table, x = snvs, by = "genome")

snvs <- snvs[order(phylum), ]
snvs <- rbind(snvs[group.ID == "Planktophila"], snvs[group.ID == "Nanopelagicus"], snvs[group.ID == "Other"], snvs[phylum != "Actinobacteriota"])
label.locs <- data.table("group.ID" = unique(snvs$group.ID), "label.loc" = 1:length(unique(snvs$group.ID)))
snvs <- merge(x = snvs, y = label.locs, by = "group.ID")

# ---- export ----

fwrite(x = snvs, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/abrupt_genomic_changes.csv")
