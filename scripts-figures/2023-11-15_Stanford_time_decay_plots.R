# RRR


# ---- set up ----

library(data.table)
library(lubridate)


# dist.file <- "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-21_3300043464_group3_bin69_all_SNV_euclidean_distance_matrix.rds" # B
# dist.file <- "data/2023-11-01_multidimensional_SNV_analysis/ME2016-07-20_3300033996_group7_bin32_all_SNV_euclidean_distance_matrix.rds" # C
dist.file <- "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-04_3300044729_group3_bin142_all_SNV_euclidean_distance_matrix.rds" # A

sample.key <- "data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv"
tax.file <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"
threads <- 1

# ---- functions ----

fill.btwn.lines <- function(X, Y1, Y2, Color, xpd = F){
  index <- is.na(X) | is.na(Y1) | is.na(Y2)
  X <- X[!index]
  Y1 <- Y1[!index]
  Y2 <- Y2[!index]
  poly.x <- c(X, X[length(X):1])
  poly.y <- c(Y1, Y2[length(Y2):1])
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

add.magic.boxes <- function(){
  box(which = "inner", col="red", lwd = 3)
  box(which = "outer", col="blue", lwd = 3)
  box(which = "plot", col="purple", lwd = 3)
  box(which = "figure", col="orange", lwd = 3)
  # inner margins (orange to purple) are mar/plt/mai
  # outer margins (red to blue) are oma/omd/omi
  par.settings <- par(no.readonly = TRUE)
  cat("change distance btwn ORANGE and PURPLE with mar  (mar = ", par.settings$mar, ")\n")
  cat("change distance btwn BLUE and RED with oma       (oma = ", par.settings$oma, ")\n")
}

# ---- import and format ----

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

# ---- get pairwise distances for time decay ----

dist.each <- as.data.table(dist.mat, keep.rownames = "sample.1")
dist.each <- melt(data = dist.each, id.vars = "sample.1", variable.name = "sample.2", value.name = "dist")
dist.each <- merge(x = dist.each, y = sample.key[ ,.(sample, date, season)], by.x = "sample.1", by.y = "sample")
dist.each <- merge(x = dist.each, y = sample.key[ ,.(sample, date, season)], by.x = "sample.2", by.y = "sample", suffixes = c(".1",".2"))
dist.each <- dist.each[ ,.(sample.1,sample.2,date.1,date.2,season.1, season.2,dist)]

dist.each$time.sec <- int_length(interval(start = dist.each[ ,date.1], end = dist.each[ ,date.2]))

dist.each <- dist.each[time.sec > 0]

dist.each[ ,time.approx.year := time.sec / 3.154e7]
dist.each[ ,time.approx.month := time.sec / 2.628e6]

dist.each <- merge(x = dist.each, y = color.key, by.x = "season.1", by.y = "season")
dist.each <- merge(x = dist.each, y = color.key, by.x = "season.2", by.y = "season", suffixes = c(".1",".2"))

dist.each <- dist.each[order(date.1, date.2)]

# calculate moving averages for time decay

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
  moving.ave[r, ave := mean(my.vals)]
  moving.ave[r, sd := sd(my.vals)]
}

moving.ave[ ,time.approx.year := center / 3.154e7]


# ---- make real-time plot ----

# pdf(file = "figures/2023-11-15_Stanford_plots/real-time_distance_acI-C.pdf", width = 6.5, height = 3)
pdf(file = "figures/2023-11-15_Stanford_plots/real-time_distance_acI-A.pdf", width = 6.5, height = 3)

par(mar = c(1.5,3.25,.25,.25), fig = c(0,1,0,1))
plot(dist ~ time.approx.year, data = dist.first[-1], type = "n", ann = F, axes = F)
points(dist ~ time.approx.year, data = dist.first[-1], pch = 19, col = adjustcolor("navy",.5))
lines(dist ~ time.approx.year, data = dist.first[-1], pch = 19, col = adjustcolor("navy",.5))
my.lm <- lm(dist ~ time.approx.year, data = dist.first[-1])
abline(my.lm)
my.lm.summ <- summary(my.lm)
R2 <- round(my.lm.summ$adj.r.squared, digits = 2)
mtext(text = bquote(adj~R^2 == ~ .(R2)), side = 3, line = -2, adj = .05)
# mtext(text = paste("pval =", round(my.lm.summ$coefficients[2,4], digits = 4)), side = 3, line = -3, adj = .1)
box()
# axis(side = 1, at = c(min(dist.first$time.approx.year), max(dist.first$time.approx.year)), labels = FALSE, lwd = 1, lwd.ticks = 0)
axis(side = 1, at = seq(0,19,2), labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 1, at = seq(0,19,2), labels = seq(2000,2019,2), lwd = 0, line = -.5)
# axis(side = 2, at = c(min(dist.first$dist), max(dist.first$dist)), labels = FALSE, lwd = 1, lwd.ticks = 0)
axis(side = 2, labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 2, las = 2, lwd = 0, line = -.35)
mtext(text = "Distance from first timepoint", side = 2, line = 2.25)
dev.off()

# ---- make time decay plot ----

# pdf(file = "figures/2023-11-15_Stanford_plots/time_decay_acI-C.pdf", width = 6.5, height = 3)
pdf(file = "figures/2023-11-15_Stanford_plots/time_decay_acI-A.pdf", width = 6.5, height = 3)
par(mar = c(2.8,3.25,.25,.25), fig = c(0,1,0,1))
plot(dist ~ time.approx.year, data = dist.each, type = "n", ann = F, axes = F)
# fill.btwn.lines(X = moving.ave$time.approx.year, 
#                 Y1 = moving.ave$ave + moving.ave$sd, 
#                 Y2 = moving.ave$ave - moving.ave$sd, Color = adjustcolor("black",.2))
points(dist ~ time.approx.year, data = dist.each, pch = 19, col = adjustcolor("royalblue",.1), cex = .5)
lines(ave ~ time.approx.year, data = moving.ave, lwd = 3, col = "black")
box()
# axis(side = 1, at = c(min(dist.first$time.approx.year), max(dist.first$time.approx.year)), labels = FALSE, lwd = 1, lwd.ticks = 0)
axis(side = 1, at = seq(0,19,2), labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 1, at = seq(0,19,2), lwd = 0, line = -.5)
# axis(side = 2, at = c(min(dist.first$dist), max(dist.first$dist)), labels = FALSE, lwd = 1, lwd.ticks = 0)
axis(side = 2, labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 2, las = 2, lwd = 0, line = -.35)
mtext(text = "Distance between timepoints", side = 2, line = 2.25)
mtext(text = "Years between timepoints", side = 1, line = 1.65)

par(mar = c(0,0,0,0), fig = c(.725,.95,.225,.4), new = TRUE)
plot(1:5,1:5, type = "n",ann = F, axes = F)
points(x = 1.25, y = 3.75, pch = 19, col = adjustcolor("royalblue",.7), cex = .5)
text(x = 1.65, y = 3.75, labels = "Pairwise comparison", adj = 0, cex = .7)
segments(x0 = 1.15, x1 = 1.35, y0 = 2, y1 = 2, lwd = 3, col = "black")
text(x = 1.65, y = 2, labels = "6 mo moving average", adj = 0, cex = .7)
box()

dev.off()
