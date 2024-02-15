# RRR
# before I very simply said "it's seasonal" if any of the PNAS seasons were different using ANOVA
# but now that I did fft to get seasonality of SNVs, that's looking very crude
# plus, the PNAS seasons seemed good for the 16S community but not as great for individual genomes, which often span 2 seasons, and had fft oscillations with only 1 season even.
# also, before I said: if it has seasonal nuc div, then do the correlation. 
# but now with the story it makes more sense to say, out of the X% of seasonally abundant ones, which also have seasonal nuc div.
# So I will use fft to define seasonally abundant ones (this script) 
# and then re-do the correlation analysis to get diverse/clonal blooms (next script)

# ---- set up ----

library(data.table)
library(lubridate)
library(compositions) # for CLR transform

fix.date <- "character"
names(fix.date) <- "date"
genome.info <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz", colClasses = fix.date)

bin.info <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds") 
bin.info <- as.data.table(bin.info)

rel.abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
colnames(rel.abund)[3] <- "abund.perc"
rel.abund <- rel.abund[Genome != "unmapped", .(sample = Sample, genome = Genome, abund.perc)]

tax <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

output.folder <- "data/2023-12-06_abundance_seasonality_analysis"

min.dates <- 30 # must be present at least 30 times, same as with SNVs but without the min coverage requirement as well
min.years <- 10 # must be present in at least 10 years, same as SNV cutoff b/c long-term analysis and SNV seasonality were in the same script. just be consistent I guess even though 10 years feels harsh

breadth.cutoff <- .7 # for instrain breadth / breadth expected
abund.cutoff <- 0 # after applying the breadth cutoff, this essentially becomes .004-.005 percent abundance (under which it's called zero)

# ---- filter out rare genomes ----

# apply breadth cutoff in addition to coverM abundance over 0
genome.info <- genome.info[(breadth / breadth_expected) >= breadth.cutoff] 
rel.abund <- rel.abund[abund.perc > abund.cutoff]
rel.abund <- merge(x = genome.info[], y = rel.abund, by = c("genome","sample"), all.x = TRUE, all.y = FALSE)

# apply # of samples cutoffs
rel.abund[ ,date := parse_date_time(date,"ymd")]
rel.abund[ ,month := month(date)] 
sample.stats <- rel.abund[ , .(temp = .N), by = .(genome, date, year, month)]
sample.dates <- sample.stats[ , .(num.dates = length(unique(date))), by = .(genome)]
sample.years <- sample.stats[ , .(num.years = length(unique(year))), by = .(genome)]
sample.months <- sample.stats[ , .(num.months = length(unique(month))), by = .(genome)] # may use later on for nucl div
sample.stats <- merge(x = sample.dates, y = sample.years, by = "genome")
sample.stats <- merge(x = sample.stats, y = sample.months, by = "genome")

sample.stats <- sample.stats[num.dates >= min.dates & num.years >= min.years]
rel.abund <- merge(x = sample.stats, y = rel.abund, by = "genome", all = FALSE)

# ---- do CLR transform on relative abundance data ----

# make a data matrix
clr.abund <- dcast(data = rel.abund, formula = genome ~ date, value.var = "abund.perc", fill = 0, fun.aggregate = mean)
clr.abund <- as.matrix(clr.abund, rownames = TRUE)
clr.abund <- clr.abund |>
  t() |>
  clr() |>
  t() |>
  as.matrix() |>
  as.data.table(keep.rownames = "genome")

# ---- interpolate and subtract big trends before fft ----

# split out each genome into elements of a list, b/c my mind is breaking thinking about all these steps with apply
genome.vals <- list()
genome.interps <- list()
for (g in clr.abund[ ,genome]){
  my.genome <- melt(data = clr.abund[genome == g], id.vars = "genome", variable.name = "date", value.name = "clr.abund")
  my.genome[ ,date := parse_date_time(date, "ymd")]
  my.genome <- my.genome[order(date)]
  my.genome[ ,time.sec := int_length(interval(start = my.genome[1,date], end = date))]
  my.genome[ ,time.approx.year := time.sec / 3.154e7]
  
  linear.fun <- approxfun(x = my.genome[ ,time.sec], y = my.genome[ ,clr.abund])
  my.interp <- data.table("time.sec" = seq.int(from = 0, to = my.genome[nrow(my.genome), time.sec], by = 86400)) # 86400 seconds in a day
  my.interp[ ,interp.clr.abund := linear.fun(v = time.sec)]
  my.interp[ ,time.approx.year := time.sec / 3.154e7]
  
  my.genome <- list(my.genome)
  names(my.genome) <- g
  genome.vals <- c(genome.vals, my.genome)
  
  my.interp <- list(my.interp)
  names(my.interp) <- g
  genome.interps <- c(genome.interps, my.interp)
}

# detrend by subtracting a cubic fit (remove any long-term trends like gradual change or big wiggles, and center around zero)
for (g in names(genome.interps)){
  my.fit <- lm(interp.clr.abund ~ time.approx.year + I(time.approx.year^2) + I(time.approx.year^3), data = genome.interps[[g]])
  my.fit <- summary(my.fit)
  my.fit <- data.table("c" = my.fit$coefficients[1,1], "x" = my.fit$coefficients[2,1], "x2" = my.fit$coefficients[3,1], "x3" = my.fit$coefficients[4,1])
  my.fit <- my.fit$c + my.fit$x * genome.interps[[g]]$time.approx.year  + my.fit$x2 * genome.interps[[g]]$time.approx.year^2 + my.fit$x3 * genome.interps[[g]]$time.approx.year^3 
  
  genome.interps[[g]][ , detrended.abund := interp.clr.abund - my.fit]
}

# ---- do fft on the detrended linear interp of CLR abundances ----

genome.ffts <- list()
for (g in names(genome.interps)){
  my.fft <- fft(z = genome.interps[[g]]$detrended.abund, inverse = F)
  my.fft <- abs(my.fft)
  my.fft <- my.fft[1:(length(my.fft) / 2)]
  
  expected.period <- 365 # x axis units are time in approximate years, but interpolated spacing is 1 day and fft only counts the number of y values, so 1 yr period is 365
  num.points <- nrow(genome.interps[[g]])
  expected.peak.index <- num.points / expected.period + 1
  
  # # manual check
  # plot(my.fft, type = "l", xlim = c(0,100))
  # abline(v = expected.peak.index, col = "red")
  
  my.fft <- data.table("index" = 1:length(my.fft), "peak" = my.fft)
  my.fft <- my.fft[order(peak, decreasing = T)]
  my.fft[ ,"period" := num.points / (index - 1)]
  
  my.fft[ ,is.year := FALSE]
  my.fft[period >= expected.period - 30 & period <= expected.period + 30 , is.year := TRUE]
  fft.stats <- data.table("genome" = g, "expected.peak.index" = expected.peak.index, "expected.period" = expected.period,
                          "peak.1.index" = my.fft[1,index], "peak.1.period" = my.fft[1,period], "peak.1.is.year" = my.fft[1,is.year],
                          "peak.2.index" = my.fft[2,index], "peak.2.period" = my.fft[2,period], "peak.2.is.year" = my.fft[2,is.year],
                          "peak.3.index" = my.fft[3,index], "peak.3.period" = my.fft[3,period], "peak.3.is.year" = my.fft[3,is.year],
                          "peak.4.index" = my.fft[4,index], "peak.4.period" = my.fft[4,period], "peak.4.is.year" = my.fft[4,is.year],
                          "peak.5.index" = my.fft[5,index], "peak.5.period" = my.fft[5,period], "peak.5.is.year" = my.fft[5,is.year])
  
  fft.stats <- list(fft.stats)
  names(fft.stats) <- g
  genome.ffts <- c(genome.ffts, fft.stats)
}

abund.ffts <- rbindlist(l = genome.ffts)

# ---- also do fft on nucleotide diversity, while we're at it ----

# make a wide table
div.tab <- dcast(data = rel.abund, formula = genome ~ date, value.var = "nucl_diversity", fill = NA, fun.aggregate = mean)

# do linear interp OVER the NA values as well (don't know div when abund is zero)
genome.vals <- list()
genome.interps <- list()
for (g in div.tab[ ,genome]){
  my.genome <- melt(data = div.tab[genome == g], id.vars = "genome", variable.name = "date", value.name = "nucl_diversity")
  my.genome <- my.genome[!is.na(nucl_diversity)]
  my.genome[ ,date := parse_date_time(date, "ymd")]
  my.genome <- my.genome[order(date)]
  my.genome[ ,time.sec := int_length(interval(start = my.genome[1,date], end = date))]
  my.genome[ ,time.approx.year := time.sec / 3.154e7]
  
  linear.fun <- approxfun(x = my.genome[ ,time.sec], y = my.genome[ ,nucl_diversity])
  my.interp <- data.table("time.sec" = seq.int(from = 0, to = my.genome[nrow(my.genome), time.sec], by = 86400)) # 86400 seconds in a day
  my.interp[ ,interp.nucl_div := linear.fun(v = time.sec)]
  my.interp[ ,time.approx.year := time.sec / 3.154e7]
  
  my.genome <- list(my.genome)
  names(my.genome) <- g
  genome.vals <- c(genome.vals, my.genome)
  
  my.interp <- list(my.interp)
  names(my.interp) <- g
  genome.interps <- c(genome.interps, my.interp)
}

# detrend by subtracting a cubic fit (remove any long-term trends like gradual change or big wiggles, and center around zero)
for (g in names(genome.interps)){
  my.fit <- lm(interp.nucl_div ~ time.approx.year + I(time.approx.year^2) + I(time.approx.year^3), data = genome.interps[[g]])
  my.fit <- summary(my.fit)
  my.fit <- data.table("c" = my.fit$coefficients[1,1], "x" = my.fit$coefficients[2,1], "x2" = my.fit$coefficients[3,1], "x3" = my.fit$coefficients[4,1])
  my.fit <- my.fit$c + my.fit$x * genome.interps[[g]]$time.approx.year  + my.fit$x2 * genome.interps[[g]]$time.approx.year^2 + my.fit$x3 * genome.interps[[g]]$time.approx.year^3 
  
  genome.interps[[g]][ , detrended.div := interp.nucl_div - my.fit]
}

# get oscillation frequencies from fast fourier transform
genome.ffts <- list()
for (g in names(genome.interps)){
  my.fft <- fft(z = genome.interps[[g]]$detrended.div, inverse = F)
  my.fft <- abs(my.fft)
  my.fft <- my.fft[1:(length(my.fft) / 2)]
  
  expected.period <- 365 # x axis units are time in approximate years, but interpolated spacing is 1 day and fft only counts the number of y values, so 1 yr period is 365
  num.points <- nrow(genome.interps[[g]])
  expected.peak.index <- num.points / expected.period + 1
  
  my.fft <- data.table("index" = 1:length(my.fft), "peak" = my.fft)
  my.fft <- my.fft[order(peak, decreasing = T)]
  my.fft[ ,"period" := num.points / (index - 1)]
  
  my.fft[ ,is.year := FALSE]
  my.fft[period >= expected.period - 30 & period <= expected.period + 30 ,is.year := TRUE]
  fft.stats <- data.table("genome" = g, "expected.peak.index" = expected.peak.index, "expected.period" = expected.period,
                          "peak.1.index" = my.fft[1,index], "peak.1.period" = my.fft[1,period], "peak.1.is.year" = my.fft[1,is.year],
                          "peak.2.index" = my.fft[2,index], "peak.2.period" = my.fft[2,period], "peak.2.is.year" = my.fft[2,is.year],
                          "peak.3.index" = my.fft[3,index], "peak.3.period" = my.fft[3,period], "peak.3.is.year" = my.fft[3,is.year],
                          "peak.4.index" = my.fft[4,index], "peak.4.period" = my.fft[4,period], "peak.4.is.year" = my.fft[4,is.year],
                          "peak.5.index" = my.fft[5,index], "peak.5.period" = my.fft[5,period], "peak.5.is.year" = my.fft[5,is.year])
  
  fft.stats <- list(fft.stats)
  names(fft.stats) <- g
  genome.ffts <- c(genome.ffts, fft.stats)
}

nuc.div.ffts <- rbindlist(l = genome.ffts)

# ---- combine and format fft stats ----

# say it's seasonal if any of the top 5 peaks are seasonal
abund.ffts[ ,`:=`(abund.fft.peak = as.numeric(NA), abund.fft.period = as.numeric(NA), abund.is.seasonal = FALSE)]
abund.ffts[peak.5.is.year == TRUE, `:=`(abund.fft.peak = 5, abund.fft.period = peak.5.period, abund.is.seasonal = TRUE)]
abund.ffts[peak.4.is.year == TRUE, `:=`(abund.fft.peak = 4, abund.fft.period = peak.4.period, abund.is.seasonal = TRUE)]
abund.ffts[peak.3.is.year == TRUE, `:=`(abund.fft.peak = 3, abund.fft.period = peak.3.period, abund.is.seasonal = TRUE)]
abund.ffts[peak.2.is.year == TRUE, `:=`(abund.fft.peak = 2, abund.fft.period = peak.2.period, abund.is.seasonal = TRUE)]
abund.ffts[peak.1.is.year == TRUE, `:=`(abund.fft.peak = 1, abund.fft.period = peak.1.period, abund.is.seasonal = TRUE)]

nuc.div.ffts[ ,`:=`(nuc.div.fft.peak = as.numeric(NA), nuc.div.fft.period = as.numeric(NA), nuc.div.is.seasonal = FALSE)]
nuc.div.ffts[peak.5.is.year == TRUE, `:=`(nuc.div.fft.peak = 5, nuc.div.fft.period = peak.5.period, nuc.div.is.seasonal = TRUE)]
nuc.div.ffts[peak.4.is.year == TRUE, `:=`(nuc.div.fft.peak = 4, nuc.div.fft.period = peak.4.period, nuc.div.is.seasonal = TRUE)]
nuc.div.ffts[peak.3.is.year == TRUE, `:=`(nuc.div.fft.peak = 3, nuc.div.fft.period = peak.3.period, nuc.div.is.seasonal = TRUE)]
nuc.div.ffts[peak.2.is.year == TRUE, `:=`(nuc.div.fft.peak = 2, nuc.div.fft.period = peak.2.period, nuc.div.is.seasonal = TRUE)]
nuc.div.ffts[peak.1.is.year == TRUE, `:=`(nuc.div.fft.peak = 1, nuc.div.fft.period = peak.1.period, nuc.div.is.seasonal = TRUE)]

genome.ffts <- merge(abund.ffts[ ,.(genome, abund.fft.peak, abund.fft.period, abund.is.seasonal)], y = nuc.div.ffts[ ,.(genome, nuc.div.fft.peak, nuc.div.fft.period, nuc.div.is.seasonal)], by = "genome", all = TRUE)
genome.ffts <- merge(x = sample.stats, y = genome.ffts, by = "genome", all.x = FALSE, all.y = TRUE)

# ---- combine clr.abunds with rel.abunds to save the data ----

clr.abund <- melt(data = clr.abund, id.vars = "genome", variable.name = "date", value.name = "clr.abund")

rel.abund <- rel.abund[ , .(nucl_diversity = mean(nucl_diversity), abund.perc = mean(abund.perc)), by = .(genome, date)]

rel.abund[ ,date := as.character(date)]
rel.abund <- merge(x = rel.abund, y = clr.abund, by = c("genome","date"), all = TRUE)
rel.abund[is.na(abund.perc), abund.perc := 0]

# ---- export data ----

fwrite(x = genome.ffts, file = file.path(output.folder,"genome_abundance_fft_stats.tsv.gz"), sep = "\t")
fwrite(x = rel.abund, file = file.path(output.folder,"filtered_abundance_and_nucl_div.tsv.gz"), sep = "\t")
