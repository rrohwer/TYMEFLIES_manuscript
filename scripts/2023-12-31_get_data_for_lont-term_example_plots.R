# RRR
# for fig 3
# modify from server-scripts/2023-11-01_multivariate_SNV_analysis/3b-make_time_decay_plots.R
# pull distance data from the server
# example file path: multivariate_SNV_analysis-med_cov_10/distance_matrices_all_SNVs/ME2015-07-03_3300042555_group6_bin161_all_SNV_euclidean_distance_matrix.rds

library(data.table)
library(lubridate)

# repeat for each file
# dist.obj <- readRDS(file = "data/2023-12-29_genome_examples_data_for_paper_figures/ME2015-07-03_3300042555_group6_bin161_all_SNV_euclidean_distance_matrix.rds") # acI-A disturbance/resilience
# my.genome <- "ME2015-07-03_3300042555_group6_bin161"
# output.file <- paste("figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/snv_dist_table-disturbance",my.genome,".csv", sep = "-")

# dist.obj <- readRDS(file = "data/2023-12-29_genome_examples_data_for_paper_figures/ME2011-09-21_3300043464_group3_bin69_all_SNV_euclidean_distance_matrix.rds") # acI-B step change
# my.genome <- "ME2011-09-21_3300043464_group3_bin69"
# output.file <-  paste("figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/snv_dist_table-step",my.genome,".csv", sep = "-")

# dist.obj <- readRDS(file = "data/2023-12-29_genome_examples_data_for_paper_figures/ME2005-06-22_3300042363_group2_bin84_all_SNV_euclidean_distance_matrix.rds") # AcAMD-5 gradual change
# my.genome <- "ME2005-06-22_3300042363_group2_bin84"
# output.file <-  paste("figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/snv_dist_table-gradual",my.genome,".csv", sep = "-")

sample.key <- fread(file = "data/2023-12-08_combined_genome_info/sample_key.tsv", colClasses = c("date" = "character"))

# ---- format into tables ----

sample.key[ ,date := parse_date_time(date, "ymd")]

dist.mat <- as.matrix(dist.obj)

# convert to long-format table
dist.each <- as.data.table(dist.mat, keep.rownames = "sample.1")
dist.each <- melt(data = dist.each, id.vars = "sample.1", variable.name = "sample.2", value.name = "dist")
dist.each <- merge(x = dist.each, y = sample.key[ ,.(sample, date, season)], by.x = "sample.1", by.y = "sample")
dist.each <- merge(x = dist.each, y = sample.key[ ,.(sample, date, season)], by.x = "sample.2", by.y = "sample", suffixes = c(".1",".2"))
dist.each <- dist.each[ ,.(sample.1,sample.2,date.1,date.2,season.1, season.2,dist)]

# get time lag between pairwise comparisons
dist.each$time.sec <- int_length(interval(start = dist.each[ ,date.1], end = dist.each[ ,date.2]))
dist.each <- dist.each[time.sec > 0]
dist.each[ ,time.approx.year := time.sec / 3.154e7]
dist.each[ ,time.approx.month := time.sec / 2.628e6]
dist.each <- dist.each[order(date.1, date.2)]

# subset to just distances from the first observation
dist.first <- dist.each[date.1 == dist.each$date.1[1]]

# ---- export file ----
dist.first[ ,genome := my.genome]
dist.first[ ,date.1 := as.character(date.1)]
dist.first[ ,date.2 := as.character(date.2)]

fwrite(x = dist.first, file = output.file)