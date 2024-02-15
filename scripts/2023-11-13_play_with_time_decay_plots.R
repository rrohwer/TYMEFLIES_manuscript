# RRR


# ---- set up ----

library(data.table)
library(lubridate)

userprefs <- commandArgs(trailingOnly = TRUE)
dist.file <- userprefs[1]
sample.key <- userprefs[2]
output.time.decay.plot <- userprefs[3]
output.directional.change.plot <- userprefs[4]
output.directional.change.stats <- userprefs[5]
threads <- as.numeric(userprefs[6])

# local paths testing
dist.file <- "data/2023-11-01_multidimensional_SNV_analysis/ME2016-07-20_3300033996_group7_bin32_all_SNV_euclidean_distance_matrix.rds"
dist.years.file <- "data/2023-11-01_multidimensional_SNV_analysis/ME2016-07-20_3300033996_group7_bin32_all_SNV_pairwise_year_stats.tsv.gz"
sample.key <- "data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv"
tax.file <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"
output.time.decay.plot.folder <- "data/2023-11-01_multidimensional_SNV_analysis/time_decay_example_plots/"
output.directional.change.plot.folder <- "data/2023-11-01_multidimensional_SNV_analysis/time_decay_example_plots/"
output.directional.change.stats <- "data/2023-11-01_multidimensional_SNV_analysis/time_decay_example_plots/ME2016-07-20_3300033996_group7_bin32_all_SNV_time_decay_stats.tsv"
threads <- 1

# ---- import and format ----

# not all genomes HAVE a dist object, they had to have enough coverage. skip if missing the file:
if (!file.exists(dist.file)){
  quit(save = "no", status = 0)
}

my.genome <- sub(pattern = "^.*ME", replacement = "", x = dist.file)
my.genome <- sub(pattern = "_all_SNV_euclidean_distance_matrix.rds", replacement = "", x = my.genome)
my.genome <- sub(pattern = "_nonsynonymous_SNV_euclidean_distance_matrix.rds", replacement = "", x = my.genome)
my.genome <- paste0("ME",my.genome)
cat("Processing", my.genome,"\n")

dist.obj <- readRDS(dist.file)

sample.key <- fread(file = sample.key, colClasses = "character")

# make sample.key match the order of samples in dist.obj
dist.order <- data.table("sample" = attributes(dist.obj)$Labels)
dist.order[ ,dist.obj.order := 1:nrow(dist.order)]
sample.key <- merge(x = dist.order, y = sample.key, by = "sample", all.x = TRUE, all.y = FALSE)
sample.key <- sample.key[order(dist.obj.order)]

# pull out genome's taxonomy
tax <- readRDS(file = tax.file)
tax <- tax[ ,c("bin.full.name","completeness","contamination","num.in.cluster","domain","phylum","class","order","family","genus","species")]
colnames(tax)[1] <- "genome"
tax <- as.data.table(tax)
tax <- tax[genome == my.genome]
tax.label <- paste(tax[ ,.(phylum,class,order,family,genus,species)], collapse = "_")

dist.years <- fread(file = dist.years.file, nThread = threads)

dist.first.year <- dist.years[year.1 == 2000]

# ---- calc distance from each sample to all first sample ----

dist.mat <- as.matrix(dist.obj)

dist.first <- data.table("sample" = names(dist.mat[1, ]), "dist" = dist.mat[1, ])

dist.first <- merge(x = dist.first, y = sample.key, by = "sample")

dist.first[ ,date := parse_date_time(x = date, orders = "ymd")]
dist.first$time.sec <- int_length(interval(start = dist.first[1,date], end = dist.first$date))
dist.first[ ,time.approx.year := time.sec / 3.154e7]

color.key <- data.table("season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), "color.season" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))
dist.first <- merge(x = dist.first, y = color.key, by = "season")

dist.first <- dist.first[order(date)]

plot(dist ~ time.approx.year, data = dist.first[-1], type = "n")
points(dist ~ time.approx.year, data = dist.first[-1], pch = 19, col = adjustcolor("navy",.5))
lines(dist ~ time.approx.year, data = dist.first[-1], pch = 19, col = adjustcolor("navy",.5))
my.lm <- lm(dist ~ time.approx.year, data = dist.first[-1])
abline(my.lm)
my.lm.summ <- summary(my.lm)
mtext(text = paste("adj R2 =", round(my.lm.summ$adj.r.squared, digits = 2)), side = 3, line = -2, adj = .1)
mtext(text = paste("pval =", round(my.lm.summ$coefficients[2,4], digits = 4)), side = 3, line = -3, adj = .1)
mtext(text = tax.label, side = 3, line = 1, cex = .7, adj = 0)
mtext(text = my.genome, side = 3, line = 2.5, adj = 0)

plot(dist ~ time.approx.year, data = dist.first[-1], type = "n")
points(dist ~ time.approx.year, data = dist.first[-1], pch = 19, col = adjustcolor(dist.first$color.season,.5))
my.lm <- lm(dist ~ time.sec, data = dist.first[-1])
abline(my.lm)
my.lm.summ <- summary(my.lm)
mtext(text = paste("adj R2 =", round(my.lm.summ$adj.r.squared, digits = 2)), side = 3, line = -2, adj = .1)
mtext(text = paste("pval =", round(my.lm.summ$coefficients[2,4], digits = 4)), side = 3, line = -3, adj = .1)
mtext(text = tax.label, side = 3, line = 1, cex = .7, adj = 0)
mtext(text = my.genome, side = 3, line = 2.5, adj = 0)

one.season <- dist.first[-1]
one.season <- one.season[season == "Late.Summer"]
plot(dist ~ time.approx.year, data = one.season, type = "n")
points(dist ~ time.approx.year, data = one.season, pch = 19, col = adjustcolor(one.season$color.season,.5))
lines(dist ~ time.approx.year, data = one.season, pch = 19, col = adjustcolor(one.season$color.season,.5))
my.lm <- lm(dist ~ time.sec, data = one.season)
abline(my.lm)
my.lm.summ <- summary(my.lm)
mtext(text = paste("adj R2 =", round(my.lm.summ$adj.r.squared, digits = 2)), side = 3, line = -2, adj = .1)
mtext(text = paste("pval =", round(my.lm.summ$coefficients[2,4], digits = 4)), side = 3, line = -3, adj = .1)
mtext(text = tax.label, side = 3, line = 1, cex = .7, adj = 0)
mtext(text = my.genome, side = 3, line = 2.5, adj = 0)

one.season <- dist.first[-1]
one.season <- one.season[season == "Fall"]
plot(dist ~ time.approx.year, data = one.season, type = "n")
points(dist ~ time.approx.year, data = one.season, pch = 19, col = adjustcolor(one.season$color.season,.5))
lines(dist ~ time.approx.year, data = one.season, pch = 19, col = adjustcolor(one.season$color.season,.5))
my.lm <- lm(dist ~ time.sec, data = one.season)
abline(my.lm)
my.lm.summ <- summary(my.lm)
mtext(text = paste("adj R2 =", round(my.lm.summ$adj.r.squared, digits = 2)), side = 3, line = -2, adj = .1)
mtext(text = paste("pval =", round(my.lm.summ$coefficients[2,4], digits = 4)), side = 3, line = -3, adj = .1)
mtext(text = tax.label, side = 3, line = 1, cex = .7, adj = 0)
mtext(text = my.genome, side = 3, line = 2.5, adj = 0)

# ---- year centroids ----

plot(dist ~ year.2, data = dist.first.year, type = "n")
points(dist ~ year.2, data = dist.first.year, pch = 19, col = adjustcolor("red4",.5))
my.lm <- lm(dist ~ year.2, data = dist.first.year)
abline(my.lm)
my.lm.summ <- summary(my.lm)
mtext(text = paste("adj R2 =", round(my.lm.summ$adj.r.squared, digits = 2)), side = 3, line = -2, adj = .1)
mtext(text = paste("pval =", round(my.lm.summ$coefficients[2,4], digits = 4)), side = 3, line = -3, adj = .1)
mtext(text = tax.label, side = 3, line = 1, cex = .7, adj = 0)
mtext(text = my.genome, side = 3, line = 2.5, adj = 0)


# ---- full time decay plots ----

dist.each <- as.data.table(dist.mat, keep.rownames = "sample.1")
dist.each <- melt(data = dist.each, id.vars = "sample.1", variable.name = "sample.2", value.name = "dist")
dist.each <- merge(x = dist.each, y = sample.key[ ,.(sample, date, season)], by.x = "sample.1", by.y = "sample")
dist.each <- merge(x = dist.each, y = sample.key[ ,.(sample, date, season)], by.x = "sample.2", by.y = "sample", suffixes = c(".1",".2"))
dist.each <- dist.each[ ,.(sample.1,sample.2,date.1,date.2,season.1, season.2,dist)]

dist.each$time.sec <- int_length(interval(start = dist.each[ ,date.1], end = dist.each[ ,date.2]))

dist.each <- dist.each[time.sec > 0]

dist.each[ ,time.approx.year := time.sec / 3.154e7]
dist.each[ ,time.approx.month := time.sec / 2.628e6]
dist.each[ ,time.month.chunk := floor(time.approx.month)]
month.avgs <- dist.each[ , .(dist.ave = mean(dist), dist.sd = sd(dist)), by = .(time.month.chunk)]

month.avgs[ ,time.approx.year := (time.month.chunk + .5) * (2.628e6/1) * (1/3.154e7)]
month.avgs <- month.avgs[order(time.approx.year)]

window <- 604800 # 1 week
window <- 2.628e6 # 1 month
window <- 7.884e6 # 3 month
window <- 1.577e7 # 6 month
points.freq <- 86399.90531424 # 1 day
moving.ave <- data.table("floor" = seq(from = 0, to = max(dist.each$time.sec), by = points.freq))
moving.ave[ ,ceiling := floor + window]
moving.ave[ ,center := floor + (window / 2)]
for(r in 1:nrow(moving.ave)){
  my.vals <- dist.each[time.sec >= moving.ave[r,floor] & time.sec <= moving.ave[r,ceiling], dist]
  moving.ave[r, ave := mean(my.vals)]
  moving.ave[r, sd := mean(my.vals)]
}

moving.ave[ ,time.approx.year := center / 3.154e7]

dist.each <- merge(x = dist.each, y = color.key, by.x = "season.1", by.y = "season")
dist.each <- merge(x = dist.each, y = color.key, by.x = "season.2", by.y = "season", suffixes = c(".1",".2"))

plot(dist ~ time.approx.year, data = dist.each, type = "n")
points(dist ~ time.approx.year, data = dist.each, pch = 19, col = adjustcolor("navy",.1), cex = .5)

lines(ave ~ time.approx.year, data = moving.ave, lwd = 3)


points(dist.ave ~ time.approx.year, data = month.avgs, pch = 19, col = "black", cex = .7)
lines(dist.ave ~ time.approx.year, data = month.avgs)
lines(dist.ave + dist.sd ~ time.approx.year, data = month.avgs)
lines(dist.ave - dist.sd ~ time.approx.year, data = month.avgs)

plot(dist ~ time.approx.year, data = dist.each, type = "n")
points(dist ~ time.approx.year, data = dist.each, pch = 19, col = adjustcolor(dist.each$color.season.1,.1), cex = .5)

plot(dist ~ time.approx.year, data = dist.each, type = "n")
points(dist ~ time.approx.year, data = dist.each, pch = 19, col = adjustcolor(dist.each$color.season.2,.1), cex = .5)

plot(dist ~ time.approx.year, data = dist.each, type = "n")
points(dist ~ time.approx.year, data = dist.each, pch = 21, col = adjustcolor(dist.each$color.season.1,.9), bg = adjustcolor(dist.each$color.season.2,.9), cex = 1)

dist.each[ ,color.season.match := "green4"]
dist.each[season.1 != season.2, color.season.match := "red3"]

plot(dist ~ time.approx.year, data = dist.each, type = "n")
points(dist ~ time.approx.year, data = dist.each, pch = 19, col = adjustcolor(dist.each$color.season.match,.1), cex = .5)
