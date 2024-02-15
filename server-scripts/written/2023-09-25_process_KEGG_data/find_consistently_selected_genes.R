# RRR

# if one strain gains against another, ALL diffs btwn strains will appear selected for
# so ID only diffs that are CONSISTENTLY selected for to try to pick out the important strain differences only

# example server call:


# ---- set up ----

library(data.table, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)

userinput <- commandArgs(trailingOnly = TRUE)
genes <- userinput[1]
annotations <- userinput[2]
coverm <- userinput[3]
genome.info <- userinput[4]
my.genome <- userinput[5]
output.file.key.genes <- userinput[6]
output.file.key.genes.by.day <- userinput[7]
threads <- as.numeric(userinput[8])

cat("Running", my.genome,"\n")

# # local files paths to test script
# genes <- "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz"
# annotations <- "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.ko-pathways_genes-unique.tsv.gz"
# coverm <- "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz"
# genome.info <- "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz"
# my.genome <- "ME2011-09-21_3300043464_group3_bin69"
# output.file.key.genes <- "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.consistently_selected_genes.tsv.gz"
# output.file.key.genes.by.day <- "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.consistently_selected_genes-by_sample.tsv.gz"
# threads <- 1

# import files
genes <- fread(file = genes, nThread = threads)
annotations <- fread(file = annotations, nThread = threads)
coverm <- fread(file = coverm, nThread = threads)
genome.info <- fread(file = genome.info, nThread = threads)

# ---- define functions ----

get.key.genes <- function(consistency.table){
  selected.genes <- consistency.table[pos > 0] # shouldn't be positive if it's not present!
  
  above.q3 <- quantile(x = selected.genes[ ,pos.perc], digits = 2)
  above.q3 <- selected.genes[pos.perc > above.q3[4], ]
  above.q3$is.Q4 <- TRUE
  
  outliers <- boxplot(selected.genes[ ,pos.perc], method = "jitter", plot = F) 
  if (length(outliers$out) > 0){
    outliers <- selected.genes[pos.perc >= min(outliers$out)]
    outliers$is.outlier <- TRUE
    
    key.genes <- merge(x = above.q3, y = outliers, all = T) 
    key.genes[is.na(is.outlier), is.outlier := FALSE]
  }else{
    key.genes <- above.q3
    key.genes$is.outlier <- FALSE 
  }
  
  return(key.genes)
}

# ---- list whether genome is there ----

colnames(coverm)[3] <- "rel.abund"
coverm <- coverm[ ,1:3]
coverm <- coverm[Genome == my.genome]
coverm[ ,date := substr(Sample, start = 3, stop = 12)]
coverm[ ,date := parse_date_time(x = date, orders = "ymd")]
coverm[ ,Sample := NULL]
coverm <- coverm[ ,.(rel.abund = mean(rel.abund)), by = date] # should choose unfiltered for the one dup day with ww and pf

genome.info <- genome.info[genome == my.genome]
genome.info[ ,date := parse_date_time(x = date, orders = "ymd")]
genome.info[ ,above.breadth.cutoff := (breadth / breadth_expected) >= .5, ]

genome <- merge(x = genome.info, y = coverm, by = "date")
genome <- genome[ ,.(genome, sample, rel.abund, above.breadth.cutoff)]
genome$genome.pres <- FALSE
genome[rel.abund > .1 & above.breadth.cutoff == TRUE, genome.pres := TRUE]

# ---- list whether gene is there and whether it's selected for ----

genes$pos <- FALSE
genes[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman), pos := TRUE]

genes$pres <- FALSE
genes[coverage >= 5 & breadth >= 1, pres := TRUE]

genes[ , sample := sub(pattern = "\\.IS_gene_info\\.tsv", replacement = "", x = sample)]
genes[ , sample := sub(pattern = "\\.gz", replacement = "", x = sample)]

genes <- merge(x = genes, genome[ , .(sample, genome.pres)], by = "sample")

# ---- tally how often are they selected vs. there? ----

consistency <- genes[ ,.(gene, sample, year, season, genome.pres, pres, pos)]

consistency$group.2012 <- "pre-2012"
consistency[year >= 2012 ,group.2012 := "post-2012"]

overall <- consistency[genome.pres == TRUE & pres == TRUE, .(pres = sum(pres), pos = sum(pos)), by = .(gene)]
overall[ ,pos.perc := pos / pres * 100]

# # too granular
# by.year <- consistency[genome.pres == TRUE, .(pres = sum(pres), pos = sum(pos)), by = .(gene, year)]
# by.year[ ,pos.perc := pos / pres * 100]

by.2012 <- consistency[genome.pres == TRUE & pres == TRUE, .(pres = sum(pres), pos = sum(pos)), by = .(gene, group.2012)]
by.2012[ ,pos.perc := pos / pres * 100]

by.season <- consistency[genome.pres == TRUE & pres == TRUE, .(pres = sum(pres), pos = sum(pos)), by = .(gene, season)]
by.season[ ,pos.perc := pos / pres * 100]

# ---- pull out consistently selected out of all the occasionally selected ----

# overall 

overall.key.genes <- get.key.genes(consistency.table = overall)

if(nrow(overall.key.genes) >= 1){
  overall.all.dates <- merge(x = overall.key.genes, y = genes[ ,-c("pos","pres")], by = "gene", all.x = TRUE, all.y = FALSE)
}


# by 2012-group 

my.groups <- unique(by.2012$group.2012)

by.2012.list <- list()
for (g in my.groups){
  key.genes <- get.key.genes(consistency.table = by.2012[group.2012 == g])
  by.2012.list <- c(by.2012.list,list(key.genes))
}

by.2012.key.genes <- rbindlist(l = by.2012.list)

if(nrow(by.2012.key.genes) >= 1){
  by.2012.all.dates <- merge(x = by.2012.key.genes, y = genes[ ,-c("pos","pres")], by = "gene", all.x = TRUE, all.y = FALSE)
}

# by season group 

my.groups <- unique(by.season$season)

by.season.list <- list()
for (g in my.groups){
  key.genes <- get.key.genes(consistency.table = by.season[season == g])
  by.season.list <- c(by.season.list,list(key.genes))
}

by.season.key.genes <- rbindlist(l = by.season.list)
colnames(by.season.key.genes)[colnames(by.season.key.genes) == "season"] <- "season.selected"

if(nrow(by.season.key.genes) >= 1){
  by.season.all.dates <- merge(x = by.season.key.genes, y = genes[ ,-c("pos","pres")], by = "gene", all.x = TRUE, all.y = FALSE)
}

# ---- make one big "key genes" file per genome ----
colnames(by.season.key.genes)[colnames(by.season.key.genes) == "season.selected"] <- "consistent.in"
colnames(by.2012.key.genes)[colnames(by.2012.key.genes) == "group.2012"] <- "consistent.in"
overall.key.genes <- data.table(overall.key.genes[ ,1], "consistent.in" = "overall", overall.key.genes[ ,-1])

key.genes <- rbind(overall.key.genes, by.2012.key.genes, by.season.key.genes)
key.genes <- merge(x = key.genes, y = annotations, by = "gene", all.x = TRUE, all.y = FALSE)

key.genes[is.na(module.description), module.description := "no KO assigned"]
key.genes[is.na(pathway.description), pathway.description := "no KO assigned"]

# ---- get daily selection stats for the key genes subset ----

key.genes.by.day <- merge(x = key.genes, y = genes[ ,-c("pres","pos")], by = "gene")

# ---- save files ----

fwrite(x = key.genes, file = output.file.key.genes, sep = "\t", nThread = threads)
fwrite(x = key.genes.by.day, file = output.file.key.genes.by.day, sep = "\t", nThread = threads)
