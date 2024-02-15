# RRR

# re-doing the analysis to get diverse/clonal blooms
# before ran correlation if had seasonal diversity pattern
# now want to run correlation if has seasonal abundance pattern
# goal: of the ones that bloom, how many have more diverse vs. less diverse blooms?

# ---- set up ----

library(data.table)
library(lubridate)
# library(readxl)

genomes <- fread(file = "data/2023-12-06_abundance_seasonality_analysis/filtered_abundance_and_nucl_div.tsv.gz", colClasses = c("date" = "character"))

seasonal.stats <- fread(file = "data/2023-12-06_abundance_seasonality_analysis/genome_abundance_fft_stats.tsv.gz")

correlation.cutoff <- .35

max.lag <- 14 # day shift to either side to allow correlations at

# ---- run Pearson correlation ----

# do exact dates correlation

for (g in unique(genomes[ ,genome])){
  cor.seasonal <- stats::cor(x = genomes[genome == g, clr.abund], y = genomes[genome == g, nucl_diversity], method = "pearson", use = "pairwise.complete.obs")
  seasonal.stats[genome == g, Pearson.correlation := cor.seasonal]
}

# do linear interpolation to do lagged correlation
# note- not adding extra time pts, correlating each timept minus the lag with the same nucl diversity values.

genomes[ ,date := parse_date_time(date, "ymd")]
genomes[ ,julian := as.numeric(julian(date))] # days since 1970-01-01

for (g in unique(genomes[ ,genome])){
  my.best.cor <- 0
  my.best.lag <- as.numeric(NA)
  my.g <- genomes[genome == g]
  my.fun <- approxfun(y = my.g$clr.abund, x = my.g$julian, method = "linear", na.rm = T)
  for (my.lag in rep(max.lag:1, each = 2) * rep(c(-1,1), times = max.lag * 2)){ # apply up to a 2-week lag
    new.dates <- my.g$julian + my.lag
    new.vals <- my.fun(v = new.dates)
    if(sum(abs(new.vals[!is.na(my.g$nucl_diversity)]), na.rm = T) > 0){
      cor.offset <- stats::cor(x = new.vals, y = my.g$nucl_diversity, method = "pearson", use = "pairwise.complete.obs")
      if (abs(cor.offset) > abs(my.best.cor)){
        my.best.cor <- cor.offset
        my.best.lag <- my.lag
      }
    }
    
  }
  seasonal.stats[genome == g, `:=`(lag.days = my.best.lag, lag.cor = my.best.cor)]
}

# ---- run linear model to see if correlated ----

for (g in unique(genomes[ ,genome])){
  lm.seasonal <- lm(nucl_diversity ~ clr.abund, data = genomes[genome == g])
  lm.seasonal <- summary(lm.seasonal)
  seasonal.stats[genome == g, adj.R2 := lm.seasonal$adj.r.squared]
  seasonal.stats[genome == g, slope := lm.seasonal$coefficients[2,1]]
  seasonal.stats[genome == g, pvalue := lm.seasonal$coefficients[2,4]]
}

# # ---- decide on cutoffs for correlated or not ----
# 
# ss <- copy(seasonal.stats)
# 
# cor.cutoff <- .3 # this really matters
# dates.cutoff <- 30 # this doesn't change much if it gets more stringent
# months.cutoff <- 1 # this doesn't change much if it gets more stringent
# seasonality.cutoff <- TRUE # seems like this is necessary, cannot replicate removing non-seasonal ones from correlated categories with the other cutoffs
# 
# ss <- ss[num.dates >= dates.cutoff & num.months >= months.cutoff]
# if (seasonality.cutoff){
#   ss <- ss[nuc.div.is.seasonal == TRUE]
# }
# ss[ ,bloom.type := ""]
# ss[abund.is.seasonal ==TRUE & Pearson.correlation > cor.cutoff , bloom.type := "More diverse"]
# ss[abund.is.seasonal ==TRUE & Pearson.correlation < -cor.cutoff , bloom.type := "Less diverse"]
# ss[abund.is.seasonal ==TRUE & Pearson.correlation >= -cor.cutoff & Pearson.correlation <= cor.cutoff , bloom.type := "No change"]
# 
# 
# data.summary <- data.table("No.change" = c(nrow(ss[bloom.type == "No change" & nuc.div.is.seasonal == FALSE]),
#                                            nrow(ss[bloom.type == "No change" & nuc.div.is.seasonal == TRUE])),
#                            "More.diverse" = c(nrow(ss[bloom.type == "More diverse" & nuc.div.is.seasonal == FALSE]),
#                                               nrow(ss[bloom.type == "More diverse" & nuc.div.is.seasonal == TRUE])),
#                            "Less.diverse" = c(nrow(ss[bloom.type == "Less diverse" & nuc.div.is.seasonal == FALSE]),
#                                               nrow(ss[bloom.type == "Less diverse" & nuc.div.is.seasonal == TRUE])))
# data.summary <- as.matrix(data.summary)
# row.names(data.summary) <- c("nuc.div.not.seasonal", "nuc.div.seasonal")
# 
# barplot(height = data.summary, legend.text = TRUE)
# 
# # use lagged correlations
# 
# ss <- copy(seasonal.stats)
# 
# ss <- ss[num.dates >= dates.cutoff & num.months >= months.cutoff]
# if (seasonality.cutoff){
#   ss <- ss[nuc.div.is.seasonal == TRUE]
# }
# ss[ ,bloom.type := ""]
# ss[abund.is.seasonal ==TRUE & lag.cor > cor.cutoff , bloom.type := "More diverse"]
# ss[abund.is.seasonal ==TRUE & lag.cor < -cor.cutoff , bloom.type := "Less diverse"]
# ss[abund.is.seasonal ==TRUE & lag.cor >= -cor.cutoff & lag.cor <= cor.cutoff , bloom.type := "No change"]
# 
# 
# data.summary <- data.table("No.change" = c(nrow(ss[bloom.type == "No change" & nuc.div.is.seasonal == FALSE]),
#                                            nrow(ss[bloom.type == "No change" & nuc.div.is.seasonal == TRUE])),
#                            "More.diverse" = c(nrow(ss[bloom.type == "More diverse" & nuc.div.is.seasonal == FALSE]),
#                                               nrow(ss[bloom.type == "More diverse" & nuc.div.is.seasonal == TRUE])),
#                            "Less.diverse" = c(nrow(ss[bloom.type == "Less diverse" & nuc.div.is.seasonal == FALSE]),
#                                               nrow(ss[bloom.type == "Less diverse" & nuc.div.is.seasonal == TRUE])))
# data.summary <- as.matrix(data.summary)
# row.names(data.summary) <- c("nuc.div.not.seasonal", "nuc.div.seasonal")
# 
# barplot(height = data.summary, legend.text = TRUE)
# 
# # # ---- make plots to choose reasonable correlation cutoff ----
# # 
# # seasonal.stats <- seasonal.stats[abund.is.seasonal == TRUE & nuc.div.is.seasonal == TRUE]
# # 
# # genomes <- merge(x = genomes, y = seasonal.stats, by = "genome", all = FALSE)
# # genomes <- genomes[order(date)]
# # 
# # for (g in unique(genomes[ ,genome])){
# #   pdf(file = paste0("figures/2023-12-07_choose_correlation_cutoff-again/", g,".pdf"), width = 6, height = 5)
# #   par(mfrow = c(2,1), mar = c(2,3,0,.5), oma = c(.1,.1,2,.1))
# #   plot(y = genomes[genome == g, nucl_diversity], x = genomes[genome == g, yday(date)], type = "p", col = adjustcolor("blue",.8), ann = F, cex = .5)
# #   mtext(text = "nuc div", side = 2, line = 2)
# #   plot(y = genomes[genome == g, clr.abund], x = genomes[genome == g, yday(date)], type = "p", col = adjustcolor("red",.8), ann = F, cex = .5)
# #   mtext(text = "clr abund", side = 2, line = 2)
# #   mtext(text = g, side = 3, outer = T, line = .5)
# #   dev.off()
# # }
# 
# # did about half manually- 150 genomes ----
# manual <- read_excel("figures/2023-12-07_choose_correlation_cutoff-again/manual_correlation_choices.xlsx")
# colnames(manual)[2:3] <- c("choice","notes")
# manual <- as.data.table(manual)
# 
# manual <- merge(x = manual, y = seasonal.stats, by = "genome")
# manual <- manual[1:150]
# manual <- manual[!is.na(choice)]
# 
# manual[choice == "n", color := "grey"]
# manual[choice == "m", color := "green"]
# manual[choice == "o", color := "red"]
# manual[ ,cex := num.dates / max(manual$num.dates) + 1]
# 
# par(mfrow = c(1,1), mar = c(.5,4,.5,.5), oma = c(.1,.1,.1,.1))
# plot(x = sample(x = 1:50, size = nrow(manual), replace = T), y = manual[ ,Pearson.correlation], type = "n", axes = F, ann = F, ylim = c(-1,1))
# box()
# axis(side = 2, las = 2, at = seq(-1,1,.25))
# mtext(text = "Pearson correlation", side = 2, line = 3)
# points(x = sample(x = 1:50, size = nrow(manual), replace = T), y = manual[ ,Pearson.correlation], pch = 21, bg = manual[ ,color], cex = manual[ ,cex])
# abline(h = 0)
# 
# # what about if I use the best offset?
# par(mfrow = c(1,1), mar = c(.5,4,.5,.5), oma = c(.1,.1,.1,.1))
# plot(x = sample(x = 1:50, size = nrow(manual), replace = T), y = manual[ ,lag.cor], type = "n", axes = F, ann = F, ylim = c(-1,1))
# box()
# axis(side = 2, las = 2, at = seq(-1,1,.25))
# mtext(text = "Pearson correlation", side = 2, line = 3)
# points(x = sample(x = 1:50, size = nrow(manual), replace = T), y = manual[ ,lag.cor], pch = 21, bg = manual[ ,color], cex = manual[ ,cex])
# abline(h = 0)
# 
# # now only look at the seasonal ones I'm including too
# par(mfrow = c(1,1), mar = c(.5,4,.5,.5), oma = c(.1,.1,.1,.1))
# plot(x = sample(x = 1:50, size = nrow(manual[nuc.div.is.seasonal == T & abund.is.seasonal == T]), replace = T), y = manual[nuc.div.is.seasonal == T & abund.is.seasonal == T,lag.cor], type = "n", axes = F, ann = F, ylim = c(-1,1))
# box()
# axis(side = 2, las = 2, at = seq(-1,1,.25))
# mtext(text = "Pearson correlation", side = 2, line = 3)
# points(x = sample(x = 1:50, size = nrow(manual), replace = T), y = manual[nuc.div.is.seasonal == T & abund.is.seasonal == T,lag.cor], pch = 21, bg = manual[nuc.div.is.seasonal == T & abund.is.seasonal == T,color], cex = manual[nuc.div.is.seasonal == T & abund.is.seasonal == T,cex])
# abline(h = 0)
# abline(h = .35 * c(-1,1))
# 
# # focus on when I "stop missing" them:
# 
# manual[choice == "m", color := adjustcolor("white", 0)]
# manual[choice == "o", color := adjustcolor("white",0)]
# 
# plot(x = sample(x = 1:50, size = nrow(manual), replace = T), y = manual[ ,Pearson.correlation], type = "n", axes = F, ann = F, ylim = c(-1,1))
# box()
# axis(side = 2, las = 2, at = seq(-1,1,.25))
# mtext(text = "Pearson correlation", side = 2, line = 3)
# points(x = sample(x = 1:50, size = nrow(manual), replace = T), y = manual[ ,Pearson.correlation], pch = 21, bg = manual[ ,color], col = manual[ ,color], cex = 2)
# abline(h = 0)
# # try to draw lines at an "amount of grey" cutoff"
# abline(h = .3)
# abline(h = -.3)
# 
# # now including offset and excluded non-seasonal
# plot(x = sample(x = 1:50, size = nrow(manual[nuc.div.is.seasonal == T & abund.is.seasonal == T]), replace = T), y = manual[nuc.div.is.seasonal == T & abund.is.seasonal == T, lag.cor], type = "n", axes = F, ann = F, ylim = c(-1,1))
# box()
# axis(side = 2, las = 2, at = seq(-1,1,.25))
# mtext(text = "Pearson correlation", side = 2, line = 3)
# points(x = sample(x = 1:50, size = nrow(manual[nuc.div.is.seasonal == T & abund.is.seasonal == T]), replace = T), y = manual[nuc.div.is.seasonal == T & abund.is.seasonal == T ,lag.cor], pch = 21, bg = manual[nuc.div.is.seasonal == T & abund.is.seasonal == T ,color], col = manual[nuc.div.is.seasonal == T & abund.is.seasonal == T ,color], cex = 2)
# abline(h = 0)
# # try to draw lines at an "amount of grey" cutoff"
# abline(h = .35)
# abline(h = -.35)
# 
# 
# par(mfrow = c(1,3), mar = c(4,4,2,.5), oma = c(.1,.1,.1,.1))
# hist(manual[choice == "n", Pearson.correlation], breaks = 50, xlim = c(-1,1))
# abline(v = .35, col = "blue", lwd = 2)
# abline(v = -.35, col = "blue", lwd = 2)
# hist(manual[choice == "m", Pearson.correlation], breaks = 25, xlim = c(-1,1), col = "green")
# abline(v = .35, col = "blue", lwd = 2)
# hist(manual[choice == "o", Pearson.correlation], breaks = 25, xlim = c(-1,1), col = "red")
# abline(v = -.35, col = "blue", lwd = 2)
# # so calling it a match/offset does not fall off a poor correlations, but missing the correlation is higher below ~ .35
# # so choose cutoff .35
# 
# par(mfrow = c(1,3), mar = c(4,4,2,.5), oma = c(.1,.1,.1,.1))
# hist(manual[choice == "n" & nuc.div.is.seasonal == T & abund.is.seasonal == T, lag.cor], breaks = 50, xlim = c(-1,1))
# abline(v = .35, col = "blue", lwd = 2)
# abline(v = -.35, col = "blue", lwd = 2)
# hist(manual[choice == "m" & nuc.div.is.seasonal == T & abund.is.seasonal == T, lag.cor], breaks = 25, xlim = c(-1,1), col = "green")
# abline(v = .35, col = "blue", lwd = 2)
# hist(manual[choice == "o" & nuc.div.is.seasonal == T & abund.is.seasonal == T, lag.cor], breaks = 25, xlim = c(-1,1), col = "red")
# abline(v = -.35, col = "blue", lwd = 2)
# # .35 still looks good


# ---- apply cutoffs and save data ----

seasonal.stats <- seasonal.stats[abund.is.seasonal == TRUE & nuc.div.is.seasonal == TRUE & Pearson.correlation >= .3, blooms.are := "More diverse"]
seasonal.stats <- seasonal.stats[abund.is.seasonal == TRUE & nuc.div.is.seasonal == TRUE & Pearson.correlation <= -.3, blooms.are := "Less diverse"]
seasonal.stats <- seasonal.stats[abund.is.seasonal == TRUE & nuc.div.is.seasonal == TRUE & Pearson.correlation < .3 & Pearson.correlation > -.3, blooms.are := "Not correlated with nucleotide diversity"]
seasonal.stats <- seasonal.stats[abund.is.seasonal == TRUE & nuc.div.is.seasonal == FALSE, blooms.are := "Nucleotide diversity is not seasonal"]

fwrite(x = seasonal.stats, file = "data/2023-12-06_abundance_seasonality_analysis/genome_fft_stats_and_bloom_diversity.tsv.gz", sep = "\t")