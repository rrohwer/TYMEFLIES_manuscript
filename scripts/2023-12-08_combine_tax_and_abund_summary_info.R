# RRR
# because I keep doing this at the start of scripts, just make a combined file of all the per-genome info

# ---- set up ----

library(data.table)
library(lubridate)
library(readxl)

fix.date <- "character"
names(fix.date) <- "date"
genome.info <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz", colClasses = fix.date)

bin.info <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds") 
bin.info <- as.data.table(bin.info)

rel.abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
colnames(rel.abund)[3] <- "abund.perc"
rel.abund <- rel.abund[Genome != "unmapped", .(sample = Sample, genome = Genome, abund.perc, Mean, RPKM, TPM)]

sample.key <- fread(file = "data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv", colClasses = fix.date)

bigtyme <- read_excel(path = "data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx")
bigtyme <- as.data.table(bigtyme)

output.file.genome.by.sample <- "data/2023-12-08_combined_genome_info/combined_genome_info_by_sample.tsv.gz"
output.file.genome <- "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv"
output.file.sample <- "data/2023-12-08_combined_genome_info/sample_key.tsv"

breadth.cutoff <- .7 # for instrain breadth / breadth expected

# ---- combine per-sample files ----

# below breadth cutoff, adjust abundance to zero
genome.info[ ,at.least.70perc.of.expected.breadth := (breadth / breadth_expected) >= breadth.cutoff] 

combo <- merge(x = genome.info, y = rel.abund, by = c("genome","sample"), all = TRUE)

combo[ ,adj.abund.perc := abund.perc]
combo[at.least.70perc.of.expected.breadth == FALSE, adj.abund.perc := 0]

# clean up taxonomy names
bin.info <- bin.info[winner == TRUE]
bin.info[phylum == "d__", domain := "Unclassified"]
bin.info[phylum == "p__", phylum := "Unclassified"]
bin.info[class == "c__", class := "Unclassified"]
bin.info[order == "o__", order := "Unclassified"]
bin.info[family == "f__", family := "Unclassified"]
bin.info[genus == "g__", genus := "Unclassified"]
bin.info[species == "s__", species := "Unclassified"]
bin.info[is.na(domain), domain := "Unclassified"]
bin.info[is.na(phylum), phylum := "Unclassified"]
bin.info[is.na(class), class := "Unclassified"]
bin.info[is.na(order), order := "Unclassified"]
bin.info[is.na(family), family := "Unclassified"]
bin.info[is.na(genus), genus := "Unclassified"]
bin.info[is.na(species), species := "Unclassified"]
bin.info[ ,domain := sub(pattern = "d__", replacement = "", x = domain)]
bin.info[ ,phylum := sub(pattern = "p__", replacement = "", x = phylum)]
bin.info[ ,class := sub(pattern = "c__", replacement = "", x = class)]
bin.info[ ,order := sub(pattern = "o__", replacement = "", x = order)]
bin.info[ ,family := sub(pattern = "f__", replacement = "", x = family)]
bin.info[ ,genus := sub(pattern = "g__", replacement = "", x = genus)]
bin.info[ ,species := sub(pattern = "s__", replacement = "", x = species)]

combo <- merge(x = combo, y = bin.info, by.x = c("genome"), by.y = c("bin.full.name"), all = TRUE)

# genome info misses samples at too low abundance
combo <- combo[ , -c("date","year","yday","season","invasion")]
combo <- merge(x = combo, y = sample.key, by = "sample", all = TRUE)

# flag duplicate samples for easy filtering to unique dates
# (the prefiltered one, the less standard depth, and all the generous donor samples (there is an actual sample from same date)):
flag.these.samples <- c("ME2002-07-17pf_3300034102","ME2006-08-18D6_3300036398",
                          "ME2018-11-08GD_3300034116","ME2018-11-08GD_3300042399","ME2018-11-08GD_3300042510","ME2018-11-08GD_3300042901","ME2018-11-08GD_3300042940","ME2018-11-08GD_3300046832")
combo[ , lesser.duplicate := FALSE]
combo[sample %in% flag.these.samples, lesser.duplicate := TRUE]

# ---- combine per-genome files ----

# get summary stats for each genome
sample.dates <- combo[lesser.duplicate == FALSE & adj.abund.perc > 0, .(num.dates = length(unique(date))), by = .(genome)]
sample.years <- combo[lesser.duplicate == FALSE & adj.abund.perc > 0, .(num.years = length(unique(year))), by = .(genome)]
sample.dates.10x <- combo[lesser.duplicate == FALSE & adj.abund.perc > 0 & coverage_median > 10, .(num.dates.cov10x = length(unique(date))), by = .(genome)]
sample.years.10x <- combo[lesser.duplicate == FALSE & adj.abund.perc > 0 & coverage_median > 10, .(num.years.cov10x = length(unique(year))), by = .(genome)]
temp <- copy(combo)
temp[ , month := month(parse_date_time(date, "ymd"))]
sample.months <- temp[lesser.duplicate == FALSE, .(month.abund = mean(adj.abund.perc)), by = .(genome, month)]
rm(temp)
sample.months <- sample.months[ , .(max.month = mean(month[month.abund == max(month.abund)])), by = .(genome)]

genome.pres <- merge(x = sample.dates, y = sample.years, by = "genome")
genome.pres <- merge(x = genome.pres, y = sample.dates.10x, by = "genome")
genome.pres <- merge(x = genome.pres, y = sample.years.10x, by = "genome")
genome.pres <- merge(x = genome.pres, y = sample.months, by = "genome")

genome.abunds <- combo[lesser.duplicate == FALSE, .(mean.abund = mean(adj.abund.perc, na.rm = T),
                                                    sum.abund = sum(adj.abund.perc, na.rm = T),
                                                    max.abund = max(adj.abund.perc, na.rm = T),
                                                    med.abund = median(adj.abund.perc, na.rm = T),
                                                    mean.nucl_diversity = mean(nucl_diversity, na.rm = T),
                                                    max.nucl_diversity = max(nucl_diversity, na.rm = T),
                                                    mean.SNV.count = mean(divergent_site_count, na.rm = T),
                                                    max.SNV.count = max(divergent_site_count, na.rm = T)), by = .(genome)]
genome.abunds <- merge(x = genome.pres, y = genome.abunds, by = "genome", all = TRUE)

colnames(bin.info)[colnames(bin.info) == "bin.full.name"] <- "genome"
genome.combo <- merge(x = bin.info, y = genome.abunds, by = "genome", all = TRUE)

num.scaffolds <- combo[!is.na(true_scaffolds), .(num_scaffolds = unique(true_scaffolds)), by = .(genome)]
genome.combo <- merge(x = genome.combo, y = num.scaffolds, by = "genome", all = TRUE)

genome.combo <- genome.combo[order(domain, phylum, class, order, family, genus, species)]

# ---- combine per-metagenome data ----

sample.key

bigtyme[ ,sample := paste0(TYMEFLIES.name,"_",taxon_oid)]
bigtyme[sample == "ME2018-06-12_NA", sample := NA]

sample.key <- merge(x = sample.key[ ,-c("date","year","yday")], y = bigtyme, by = "sample", all = TRUE)
sample.key[ ,date := paste(Year,Month,Day, sep = "-")]
sample.key[ ,yday := yday(parse_date_time(date,"ymd"))]

sample.key[ , lesser.duplicate := FALSE]
sample.key[sample %in% flag.these.samples, lesser.duplicate := TRUE]

# ---- export data ----

fwrite(x = combo, file = output.file.genome.by.sample, sep = "\t")
fwrite(x = genome.combo, file = output.file.genome, sep = "\t")
fwrite(x = sample.key, file = output.file.sample, sep = "\t")















