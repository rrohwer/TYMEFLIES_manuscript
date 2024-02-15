# RRR
# get a clen data table for the time decay figure
# this is fig2A
# modify from server-scripts/2023-11-01_multivariate_SNV_analysis/3b-make_time_decay_plots.R
# pull distance data from the server
# example genome choice is: ME2017-06-13_3300043469_group7_bin14
# file path: multivariate_SNV_analysis-med_cov_10/distance_matrices_all_SNVs/ME2017-06-13_3300043469_group7_bin14_all_SNV_euclidean_distance_matrix.rds

library(data.table)
library(lubridate)

dist.obj <- readRDS(file = "data/2023-12-29_genome_examples_data_for_paper_figures/ME2017-06-13_3300043469_group7_bin14_all_SNV_euclidean_distance_matrix.rds")
sample.key <- fread(file = "data/2023-12-08_combined_genome_info/sample_key.tsv", colClasses = c("date" = "character"))

# ----

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

# calculate moving average to show a smoothed line in figure
# window <- 604800 # 1 week
# window <- 2.628e6 # 1 month
# window <- 7.884e6 # 3 month
window <- 1.577e7 # 6 month
points.freq <- 86399.90531424 # 1 day
moving.ave <- data.table("floor" = seq(from = 0, to = max(dist.each$time.sec), by = points.freq))
moving.ave[ ,ceiling := floor + window]
moving.ave[ ,center := floor + (window / 2)]
for(r in 1:nrow(moving.ave)){
  my.vals <- dist.each[time.sec >= moving.ave[r,floor] & time.sec <= moving.ave[r,ceiling], dist]
  moving.ave[r, dist := mean(my.vals)]
  moving.ave[r, sd := sd(my.vals)]
  moving.ave[r, num.obs := length(my.vals)]
}
moving.ave[ ,time.approx.year := center / 3.154e7]

# ---- save data to make the plot ----
dist.each[ ,`:=`(date.1 = as.character(date.1), date.2 = as.character(date.2))]

fwrite(x = dist.each, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/time_decay_distances-ME2017-06-13_3300043469_group7_bin14.csv.gz")
fwrite(x = moving.ave, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/time_decay_6mo_moving_average-ME2017-06-13_3300043469_group7_bin14.csv.gz")

# ----



