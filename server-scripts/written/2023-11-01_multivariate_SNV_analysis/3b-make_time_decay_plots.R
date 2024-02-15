# RRR

# 1. make real-time distance from first date plot to look for directional change
# fit linear lines to check for gradual directional change
#   try both fitting exact, uneven points
#   and fitting an evenly spaced linear interpolation, since extra summer samples at the end of the time series drag the line to one side
# 2. make time-decay pairwiase distances to look for seasonality
#   take 6 month moving average to draw a smooth line overlayed on the plot
# check that main period is 1 year using a fast fourrier transform (FFT)
#   linearly interpolate averaging duplicate comparisons and fit a cubic line to subtract big trends before doing FFT
#   pull out the main period- is it close to 1 year oscillations?
# 3. save one combined plot and one one-line data file for each run

# ---- set up ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(strucchange))

userprefs <- commandArgs(trailingOnly = TRUE)
dist.file <- userprefs[1]
sample.key <- userprefs[2]
tax.file <- userprefs[3]
output.plot.folder <- userprefs[4]
output.stats.folder <- userprefs[5]
threads <- as.numeric(userprefs[6])

# # local paths troubleshooting errors:
# dist.file <- "data/2023-11-01_multidimensional_SNV_analysis/example_input_files/ME2002-11-07_3300042865_group1_bin56_all_SNV_euclidean_distance_matrix.rds"

# # files used to build classifier:
# step-change
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2003-07-11_3300042290_group2_bin162_all_SNV_euclidean_distance_matrix.rds"
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2011-09-21_3300043464_group3_bin69_all_SNV_euclidean_distance_matrix.rds"
# output.plot.folder <- "data/2023-12-04_long-term_change_classifier/test_output_files/step"
# output.stats.folder <- "data/2023-12-04_long-term_change_classifier/test_output_files/step"
# # disturbance
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2004-08-06_3300043437_group2_bin96_all_SNV_euclidean_distance_matrix.rds"
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2005-06-22_3300042363_group2_bin34_all_SNV_euclidean_distance_matrix.rds"
# output.plot.folder <- "data/2023-12-04_long-term_change_classifier/test_output_files/dist"
# output.stats.folder <- "data/2023-12-04_long-term_change_classifier/test_output_files/dist"
# gradual change
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2009-06-09_3300042893_group3_bin124_all_SNV_euclidean_distance_matrix.rds"
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2005-06-22_3300042363_group2_bin84_all_SNV_euclidean_distance_matrix.rds"
# output.plot.folder <- "data/2023-12-04_long-term_change_classifier/test_output_files/grad"
# output.stats.folder <- "data/2023-12-04_long-term_change_classifier/test_output_files/grad"
# NO CHANGE
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2018-04-25_3300042375_group8_bin150_all_SNV_euclidean_distance_matrix.rds"
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2002-11-07_3300042865_group1_bin56_all_SNV_euclidean_distance_matrix.rds"
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2017-06-24_3300042306_group7_bin66_all_SNV_euclidean_distance_matrix.rds" # this is a breakpoint, check not a disturbance though
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2001-09-10pf_3300042914_group1_bin12_all_SNV_euclidean_distance_matrix.rds"
# dist.file <- "data/2023-12-04_long-term_change_classifier/test_input_files/ME2005-06-22_3300042363_group2_bin56_all_SNV_euclidean_distance_matrix.rds"
# output.plot.folder <- "data/2023-12-04_long-term_change_classifier/test_output_files/flat"
# output.stats.folder <- "data/2023-12-04_long-term_change_classifier/test_output_files/flat"
# 
# sample.key <- "data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv"
# tax.file <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"
# threads <- 1

disturbance.duration.cutoff <- 2.628e+6 # 2.628e+6 seconds in a month

# ---- define functions ----

get.genome <- function(dist.file){
  my.genome <- sub(pattern = "^.*ME", replacement = "", x = dist.file)
  my.genome <- sub(pattern = "_all_SNV_euclidean_distance_matrix.rds", replacement = "", x = my.genome)
  my.genome <- sub(pattern = "_nonsynonymous_SNV_euclidean_distance_matrix.rds", replacement = "", x = my.genome)
  my.genome <- paste0("ME",my.genome)
  return(my.genome)
}

make.sample.key.match <- function(dist.obj, sample.key){
  dist.order <- data.table("sample" = attributes(dist.obj)$Labels)
  dist.order[ ,dist.obj.order := 1:nrow(dist.order)]
  sample.key <- merge(x = dist.order, y = sample.key, by = "sample", all.x = TRUE, all.y = FALSE)
  sample.key <- sample.key[order(dist.obj.order)]
  return(sample.key)
}

get.taxonomy.vector <- function(tax){
  tax <- tax[ ,c("bin.full.name","completeness","contamination","num.in.cluster","domain","phylum","class","order","family","genus","species")]
  colnames(tax)[1] <- "genome"
  tax <- as.data.table(tax)
  tax <- tax[genome == my.genome]
  tax.label <- paste(tax[ ,.(phylum,class,order,family,genus,species)], collapse = "_")
  return(tax.label)
}

get.dist.data.table <- function(dist.mat){
  dist.each <- as.data.table(dist.mat, keep.rownames = "sample.1")
  dist.each <- melt(data = dist.each, id.vars = "sample.1", variable.name = "sample.2", value.name = "dist")
  dist.each <- merge(x = dist.each, y = sample.key[ ,.(sample, date, season)], by.x = "sample.1", by.y = "sample")
  dist.each <- merge(x = dist.each, y = sample.key[ ,.(sample, date, season)], by.x = "sample.2", by.y = "sample", suffixes = c(".1",".2"))
  dist.each <- dist.each[ ,.(sample.1,sample.2,date.1,date.2,season.1, season.2,dist)]
  return(dist.each)
}

add.time.intervals <- function(dist.table){
  # note also removes diagonal and upper half of distance matrix (so no duplicate pairwise comparisons)
  dist.table$time.sec <- int_length(interval(start = dist.table[ ,date.1], end = dist.table[ ,date.2]))
  
  dist.table <- dist.table[time.sec > 0]
  
  dist.table[ ,time.approx.year := time.sec / 3.154e7]
  dist.table[ ,time.approx.month := time.sec / 2.628e6]
  
  dist.table <- dist.table[order(date.1, date.2)]
  
  return(dist.table)
}

get.moving.average <- function(window, points.freq, dist.table){
  # window <- 604800 # 1 week
  # window <- 2.628e6 # 1 month
  # window <- 7.884e6 # 3 month
  # window <- 1.577e7 # 6 month
  # points.freq <- 86399.90531424 # 1 day
  moving.ave <- data.table("floor" = seq(from = 0, to = max(dist.table$time.sec), by = points.freq))
  moving.ave[ ,ceiling := floor + window]
  moving.ave[ ,center := floor + (window / 2)]
  for(r in 1:nrow(moving.ave)){
    my.vals <- dist.table[time.sec >= moving.ave[r,floor] & time.sec <= moving.ave[r,ceiling], dist]
    moving.ave[r, dist := mean(my.vals)]
    moving.ave[r, sd := mean(my.vals)]
  }
  moving.ave[ ,time.approx.year := center / 3.154e7]
  return(moving.ave)
}

get.linear.interpolation <- function(dist.table, points.freq){
  linear.fun <- approxfun(x = dist.table[ ,time.sec], y = dist.table[ ,dist], method = "linear", ties = mean)
  new.x <- seq(from = min(dist.table[ ,time.sec]), to = max(dist.table[ ,time.sec]), by = points.freq)
  linear.interp.dist.table <- data.table("time.sec" = new.x, "dist" = linear.fun(v = new.x))
  linear.interp.dist.table[ ,time.approx.year := time.sec / 3.154e7]
  return(linear.interp.dist.table)
}

get.linear.model.stats <- function(my.lm){
  my.lm <- summary(my.lm)
  my.lm.stats <- data.table("adj.R2" = my.lm$adj.r.squared, "pval" = my.lm$coefficients[2,4], "slope" = my.lm$coefficients[2,1], "y.intercept" = my.lm$coefficients[1,1])
  return(my.lm.stats)
}

# for paper: I identified annual oscillations in SNV frequencies using a periodogram spectral analysis by computing the 
# absolute value of a fast Fourier transform (taken using the fft function in R).
# I took a linear interpolation of all averaged pairwise time interval distances to create evenly spaced observations,
# and performed a trend removal using a cubic fit. I considered an annual pattern present if one of the top XX peaks 
# in the periodogram spectrum occurred within a month of a 1-year interval. 

get.cubic.fit.model.stats <- function(my.fit){
  my.fit <- summary(my.fit)
  my.fit.stats <- data.table("c" = my.fit$coefficients[1,1], "x" = my.fit$coefficients[2,1], "x2" = my.fit$coefficients[3,1], "x3" = my.fit$coefficients[4,1])
  return(my.fit.stats)
}

get.fft <- function(dist.table, stats.table){
  my.y <- dist.table$dist
  my.fit <- stats.table$c + stats.table$x * dist.table$time.approx.year  + stats.table$x2 * dist.table$time.approx.year^2 + stats.table$x3 * dist.table$time.approx.year^3 
  my.y <- my.y - my.fit
  my.fft <- fft(z = my.y, inverse = F)
  my.fft <- abs(my.fft)
  my.fft <- my.fft[1:(length(my.fft) / 2)]
  return(my.fft)
}

get.period.etc.from.fft <- function(my.fft, expected.period, num.points){
  # x axis units are time in approximate years, but interpolated spacing is 1 day and fft only counts the number of y values, so 1 yr period is 365
  expected.peak.index <- num.points / expected.period + 1
  my.fft <- data.table("index" = 1:length(my.fft), "peak" = my.fft)
  my.fft <- my.fft[order(peak, decreasing = T)]
  my.fft[ ,"period" := num.points / (index - 1)]
  my.fft[ ,is.year := FALSE]
  my.fft[period >= expected.period - 30 & period <= expected.period + 30 ,is.year := TRUE]
  fft.stats <- data.table("expected.peak.index" = expected.peak.index, "expected.period" = expected.period,
                          "peak.1.index" = my.fft[1,index], "peak.1.period" = my.fft[1,period], "peak.1.is.year" = my.fft[1,is.year],
                          "peak.2.index" = my.fft[2,index], "peak.2.period" = my.fft[2,period], "peak.2.is.year" = my.fft[2,is.year],
                          "peak.3.index" = my.fft[3,index], "peak.3.period" = my.fft[3,period], "peak.3.is.year" = my.fft[3,is.year],
                          "peak.4.index" = my.fft[4,index], "peak.4.period" = my.fft[4,period], "peak.4.is.year" = my.fft[4,is.year],
                          "peak.5.index" = my.fft[5,index], "peak.5.period" = my.fft[5,period], "peak.5.is.year" = my.fft[5,is.year])
  return(fft.stats)
}

# NEW BREAKPOINT ANALYSIS

get.best.breakpoint <- function(linear.interp){
  my.fstats <<- Fstats(formula = linear.interp$dist ~ 1) # need this object for making the plot, so define into global env.
  breakpoint.index <- my.fstats$breakpoint # returns whether signif or not, this is just the peak
  
  p.val <- sctest(my.fstats) # p-value for whether the lines on either side of the peak are different
  p.val <- p.val$p.value
  
  if (p.val < 0.01){
    my.breakpoint <- linear.interp.dist.first$time.approx.year[breakpoint.index]  
  }else{
    my.breakpoint <- NA
  }
  return(my.breakpoint)
}

rule.out.if.step.too.little <- function(actual.data, break.loc){
  first.half <- actual.data[time.approx.year < break.loc]
  second.half <- actual.data[time.approx.year > break.loc]
  
  if (nrow(first.half) > 0 & nrow(second.half) > 0){ # since break is on interpolated data, could be missing actual data points
    # do t-test (Mann-Whitney) on true data points, so that also reflects uncertainty from less frequent data than the interp values
    my.p <- wilcox.test(x = first.half$dist, y = second.half$dist)
    my.p <- my.p$p.value 
    
    # in addition to being statistically different, the after mean should be meaningfully/obviously higher
    before.mean <- mean(first.half$dist)
    after.mean <- mean(second.half$dist)
    step.diff <- before.mean / after.mean # before / after must be < .75
    
  }else{
    my.p <- 1
  }
  
  if (my.p >= 0.01 | step.diff >= .75){
    break.loc <- NA 
  }
  
  return(break.loc)
}

get.step.change.stats <- function(break.loc, first.date){
  start.year <- parse_date_time(first.date, "ymd") |>
    year()
  
  if (!is.na(break.loc)){
    change.date <- date_decimal(break.loc + start.year) |>
      round_date(unit = "day") |>
      as.character()
  }else{
    change.date <- NA
  }
  
  return(data.table("breakpoint.loc" = break.loc, "breakpoint.date" = change.date))
}

find.outliers <- function(dist.first){
  outliers <- boxplot(dist.first$dist, plot = F)
  outliers <- outliers$out[outliers$out > outliers$stats[3,1]] # only outliers that are high (more different), not low (more similar to first sample)
  dist.first$outlier <- FALSE
  dist.first[dist %in% outliers, outlier := TRUE]
  
  # get starts and ends of consecutive groups of outliers
  dist.first$group.start <- FALSE
  dist.first$group.end <- FALSE
  for (r in 2:nrow(dist.first)){
    if (dist.first[r, outlier] == TRUE & dist.first[r - 1, outlier] == FALSE){
      dist.first[r, group.start := TRUE]
    }else if(dist.first[r, outlier] == FALSE & dist.first[r - 1, outlier] == TRUE){
      dist.first[r - 1, group.end := TRUE]
    }
  }
  
  # make all outlier groups have a start and an end
  if(dist.first[nrow(dist.first), outlier] == TRUE & dist.first[nrow(dist.first), group.start == TRUE]){ # if ends with outlier of group length 1, just remove it now
    dist.first[nrow(dist.first), group.start := FALSE]
  }else if(dist.first[nrow(dist.first), outlier] == TRUE){ # if outlier group continues to the end, call the last point the end of that group
    dist.first[nrow(dist.first), group.end := TRUE]
  }
  
  # find paired group start and stop indexes
  dist.first[ ,group.start := cumsum(group.start)]
  dist.first[duplicated(group.start), group.start := 0]
  dist.first[ ,group.end := cumsum(group.end)]
  dist.first[duplicated(group.end), group.end := 0]
  
  return(dist.first)
}

get.table.of.disturbances <- function(dist.first, disturbance.duration.cutoff){
  # make a table of all the disturbances that are long enough
  num.outlier.groups <- max(dist.first$group.start)
  if(num.outlier.groups > 0){
    disturbances.list <- list()
    for (n in 1:(num.outlier.groups)){
      my.group <- dist.first[which(group.start == n):which(group.end == n)]
      my.duration <- my.group[nrow(my.group), time.sec] - my.group[1, time.sec]
      my.recovery <- dist.first[nrow(dist.first), time.sec] - my.group[nrow(my.group), time.sec]
      if (my.duration > disturbance.duration.cutoff & my.recovery > 2.628e+6){ # have a least a month recovery before the time series end
        my.disturbance <- data.table("start.loc" = my.group[1, time.approx.year],
                                     "end.loc" = my.group[nrow(my.group), time.approx.year],
                                     "start.date" = my.group[1, date.2],
                                     "end.date" =  my.group[nrow(my.group), date.2],
                                     "duration.days" = my.duration / 86400,
                                     "max.distance" = max(my.group[ ,dist]))
        disturbances.list <- c(disturbances.list, list(my.disturbance))
      }
    }
    disturbances.table <- rbindlist(disturbances.list)
    if(nrow(disturbances.table) < 1){
      disturbances.table <- data.table("start.loc" = NA, "end.loc" = NA, "start.date" = NA, "end.date" =  NA, "duration.days" = NA, "max.distance" = NA)
    }
  }else{
    disturbances.table <- data.table("start.loc" = NA, "end.loc" = NA, "start.date" = NA, "end.date" =  NA, "duration.days" = NA, "max.distance" = NA)
  }
  
  return(disturbances.table)
}

get.outlier.stats <- function(disturbances.table){
  my.max <- which(disturbances.table[ ,max.distance] == max(disturbances.table[ ,max.distance]))
  colnames(disturbances.table) <- paste0("max.disturbance.", colnames(disturbances.table))
  if (!is.na(disturbances.table[1,1])){
    my.stats <- data.table("total.disturbances" = nrow(disturbances.table),
                           disturbances.table[my.max, -c("max.disturbance.max.distance")])
  }else{
    my.stats <- data.table("total.disturbances" = 0,
                           disturbances.table[1, -c("max.disturbance.max.distance")])
  }
  
  return(my.stats)
}

# ---- import and format ----

my.genome <- get.genome(dist.file = dist.file)

# not all genomes HAVE a dist object, they had to have enough coverage. skip if missing the file:
if (!file.exists(dist.file)){
  quit(save = "no", status = 0)
}

dist.obj <- readRDS(dist.file)

sample.key <- fread(file = sample.key, colClasses = "character")
sample.key <- make.sample.key.match(dist.obj = dist.obj, sample.key = sample.key)

# skip genomes that are present on fewer than 30 dates 
# and present in fewer than 10 years
if (length(unique(sample.key$date)) < 30 | length(unique(sample.key$year)) < 10 ){
  quit(save = "no", status = 0)
}

cat("Processing", my.genome,"\n")

tax <- readRDS(file = tax.file)
tax.label <- get.taxonomy.vector(tax = tax)

# ---- pull distances from distance object ----

dist.mat <- as.matrix(dist.obj)

dist.each <- get.dist.data.table(dist.mat = dist.mat)

dist.each <- add.time.intervals(dist.table = dist.each)

dist.first <- dist.each[sample.1 == dist.each[1, sample.1]]

# ---- calc smoothed and interpolated time points at each day ----

moving.ave.dist.each <- get.moving.average(window = 1.577e7, points.freq = 86399.90531424, dist.table = dist.each) # 6 mo, for time decay visualization

linear.interp.dist.each <- get.linear.interpolation(dist.table = dist.each, points.freq = 86399.90531424) # daily, for FFT calculation

dist.first.trim <- copy(dist.first)
dist.first.trim <- dist.first.trim[time.sec > 2.628e+6] # minus points closer than a month to the 1st point since they skew way down

linear.interp.dist.first <- get.linear.interpolation(dist.table = dist.first.trim, points.freq = 86399.90531424) # daily, for long-term patterns, but minus points closer than a month to the 1st point since they skew way down

# ---- get directional change ----

# simple fit to all data points
lm.dist.first <- lm(dist ~ time.approx.year, data = dist.first)
lm.dist.first.stats <- get.linear.model.stats(my.lm = lm.dist.first)

# fit to interp data points to account for uneven sampling, minus first 5 data points that tend to skew low because very close to starting point
lm.dist.first.interp <- lm(dist ~ time.approx.year, data = linear.interp.dist.first)
lm.dist.first.interp.stats <- get.linear.model.stats(my.lm = lm.dist.first.interp)

# ---- get big scale wiggly and directional change to subtract before fft ----

fit.dist.each.interp <- lm(dist ~ time.approx.year + I(time.approx.year^2) + I(time.approx.year^3), data = linear.interp.dist.each) 
fit.dist.each.interp.stats <- get.cubic.fit.model.stats(my.fit = fit.dist.each.interp)

# ---- do fourrier transform ----

# dist from first is too rough, using a moving average introduces artifacts, use the linear interp of the time decay
fft.dist.each <- get.fft(dist.table = linear.interp.dist.each, stats.table = fit.dist.each.interp.stats)
fft.dist.each.stats <- get.period.etc.from.fft(my.fft = fft.dist.each, expected.period = 365, num.points = nrow(linear.interp.dist.each))

# ---- do breakpoint analysis for abrupt step change ----

breakpoint.loc <- get.best.breakpoint(linear.interp = linear.interp.dist.first) # speed up by subsampling the linear interp to less frequent values
breakpoint.loc <- rule.out.if.step.too.little(actual.data = dist.first.trim, break.loc = breakpoint.loc) # remove first data points as did with linear interp too

step.change.stats <- get.step.change.stats(break.loc = breakpoint.loc, first.date = dist.first$date.1[1])

# ---- do outlier analysis for disturbance/resilience pattern ----

dist.first.trim <- find.outliers(dist.first = dist.first.trim)
disturbances.table <- get.table.of.disturbances(dist.first = dist.first.trim, disturbance.duration.cutoff = disturbance.duration.cutoff)
disturbances.stats <- get.outlier.stats(disturbances.table = disturbances.table)

# ---- CLASSIFY and save stats ----

temp <- copy(lm.dist.first.interp.stats)
colnames(temp) <- paste("interp",colnames(temp), sep = ".") 
genome.stats <- data.table("genome" = my.genome, lm.dist.first.stats, temp, fft.dist.each.stats, step.change.stats, disturbances.stats)

genome.stats[ ,Classified.Seasonal := FALSE]
genome.stats[peak.1.is.year == TRUE | peak.2.is.year == TRUE | peak.3.is.year == TRUE | peak.4.is.year == TRUE | peak.5.is.year == TRUE, Classified.Seasonal := TRUE]

if (genome.stats$interp.adj.R2 >= .55){
  genome.stats$Classified.LT.Change <- "gradual"
}else if(!is.na(genome.stats$breakpoint.loc)){
  genome.stats$Classified.LT.Change <- "step"
}else if(genome.stats$total.disturbances > 0){
  genome.stats$Classified.LT.Change <- "disturbance"
}else{
  genome.stats$Classified.LT.Change <- "none"
}

fwrite(x = genome.stats, file = file.path(output.stats.folder,paste0(my.genome,"_SNV_freq_time_decay_stats.tsv")))

# ---- make plot ----

pdf(file = file.path(output.plot.folder,paste0(tax.label," ",my.genome,".pdf")), width = 10, height = 6)

par(mar = c(1.5,3.25,1.55,.25), fig = c(0,.75,.49,.98))
plot(dist ~ time.approx.year, data = dist.first, type = "n", ann = F, axes = F)
points(dist ~ time.approx.year, data = dist.first, pch = 19, col = adjustcolor("navy",.5))
lines(dist ~ time.approx.year, data = dist.first, col = adjustcolor("navy",.3))
abline(lm.dist.first, col = "black")
abline(lm.dist.first.interp, col = "red3")
R2 <- round(lm.dist.first.stats$adj.R2, digits = 2)
slope <- round(lm.dist.first.stats$slope, digits = 2)
mtext(text = bquote(adj~R^2 == ~ .(R2)), side = 3, line = -1.25, adj = .05, col = "black", cex = .7)
mtext(text = bquote(slope == ~ .(slope)), side = 3, line = -2, adj = .05, col = "black", cex = .7)
R2 <- round(lm.dist.first.interp.stats$adj.R2, digits = 2)
slope <- round(lm.dist.first.interp.stats$slope, digits = 2)
mtext(text = bquote(adj~R^2 == ~ .(R2)), side = 3, line = -1.25, adj = .25, col = "red3", cex = .7)
mtext(text = bquote(slope == ~ .(slope)), side = 3, line = -2, adj = .25, col = "red3", cex = .7)
axis(side = 1, labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 1, labels = TRUE, lwd = 0, line = -.5)
axis(side = 2, labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 2, las = 2, lwd = 0, line = -.35)
mtext(text = "Distance from first timepoint", side = 2, line = 2.25)
mtext(text = tax.label, side = 3, line = 1, cex = .7, adj = 0)
mtext(text = my.genome, side = 3, line = .25, cex = .7, adj = 0)
mtext(text = genome.stats$Classified.LT.Change, side = 3, line = -1.5, cex = 2, adj = .9, col = "red3", outer = TRUE)
box()
if (!is.na(disturbances.stats$max.disturbance.start.loc)){
  axis(side = 1, at = c(disturbances.stats$max.disturbance.start.loc, disturbances.stats$max.disturbance.end.loc),col = "orange2", lwd = 8, lwd.ticks = 0, labels = F)  
}
if (!is.na(step.change.stats$breakpoint.loc)){
  abline(v = step.change.stats$breakpoint.loc, col = adjustcolor("purple2",.5), lwd = 3) 
}

par(mar = c(2.8,3.25,.25,.25), fig = c(0,.75,0,.49), new = T)
plot(dist ~ time.approx.year, data = dist.each, type = "n", ann = F, axes = F)
points(dist ~ time.approx.year, data = dist.each, pch = 19, col = adjustcolor("royalblue",.1), cex = .5)
lines(dist ~ time.approx.year, data = moving.ave.dist.each, lwd = 3, col = "black")
box()
axis(side = 1, labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 1, lwd = 0, line = -.5)
axis(side = 2, labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 2, las = 2, lwd = 0, line = -.35)
mtext(text = "Distance between timepoints", side = 2, line = 2.25)
mtext(text = "Years Apart", side = 1, line = 1.65)

par(mar = c(2.8,2,.25,.25), fig = c(.75,.875,0,.49), new = T)
x.lim <- c(0, max(fft.dist.each.stats$expected.peak.index, fft.dist.each.stats$actual.peak.index) + max(fft.dist.each.stats$expected.peak.index, fft.dist.each.stats$actual.peak.index) / 2)
plot(fft.dist.each, type = "n", ann = F, axes = F)
box()
abline(v = x.lim, col = "magenta")
lines(fft.dist.each)
mtext(text = "Spectrum", side = 2, line = .5)

par(mar = c(2.8,0,.25,.25), fig = c(.875,1,0,.49), new = T)
plot(fft.dist.each, xlim = x.lim, type = "n", ann = F, axes = F)
abline(v = fft.dist.each.stats$expected.peak.index, col = "green4", lwd = 1)
# abline(v = fft.dist.each.stats$actual.peak.index, col = "royalblue", lwd = 1)
lines(fft.dist.each)
box(col = "magenta")
mtext(text = c("365 days", paste(round(fft.dist.each.stats$peak.1.period, digits = 1),"days")), side = 3, line = c(-2,-3), at = 0, adj = 0, col = c("green4", "royalblue"), cex = .7)
mtext(text = "Frequency", side = 1, line = 1, at = -.5)
mtext(text = "zoom", side = 1, line = -1, col = "magenta", adj = .2)

par(mar = c(1.5,0,1.55,.25), fig = c(.75,.8,.49,.98), new = T)
boxplot(x = dist.first$dist, col = adjustcolor("navy",.5), border = "navy", lty = 1, axes = F, ann = F, outpch = 19, outcol =  adjustcolor("navy",.5), xlim = c(0,2), at = .5)
box()
if(disturbances.stats$total.disturbances == 1){
  mtext(text = paste(disturbances.stats$total.disturbances, "disturbance"), side = 4, line = -1, col = "orange2")
}else{
  mtext(text = paste(disturbances.stats$total.disturbances, "disturbances"), side = 4, line = -1, col = "orange2")
}

par(mar = c(1.5,1,1.55,.25), fig = c(.8,1,.49,.98), new = T)
plot(my.fstats, ann = F, axes = F)
abline(v = my.fstats$breakpoint / my.fstats$nobs, col = adjustcolor("purple2",.5))
box()
mtext(text = "Breakpoint Analysis", side = 1, col = "purple2")
if (!is.na(genome.stats$breakpoint.loc)){
  mtext(text = genome.stats$breakpoint.date, side = 3, adj = .2, line = -2, col = "purple2")
}else{
  mtext(text = "No breakpoint", side = 2, adj = .9, line = -1, col = "purple2")
}


dev.off()

# ---- end ----