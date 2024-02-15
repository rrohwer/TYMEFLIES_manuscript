# RRR

# ---- set up ----

library(data.table)

seas <- fread(file = "data/2023-12-06_abundance_seasonality_analysis/genome_fft_stats_and_bloom_diversity.tsv.gz")
tax <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv")
snv <- fread(file = "data/2023-12-04_long-term_change_classifier/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")

# ---- little stats ----

nrow(seas) / nrow(tax) * 100 # 51.62872 % of genomes appeared at least 30 times and in at least 10 years

nrow(seas[abund.is.seasonal == TRUE]) / nrow(seas) * 100 # of those, 71.981 % had a seasonal abundance pattern

nrow(seas[nuc.div.is.seasonal == TRUE]) / nrow(seas) * 100 # and 32.97151 % had a seasonal nucleotide diversity pattern

nrow(snv) / nrow(tax) * 100 # 9.211909 % of genomes appeared at least 30 times at median coverage >10 and in at least 10 years

nrow(snv[Classified.Seasonal == TRUE]) / nrow(snv) * 100 # of those, 79.84791 % had seasonal strain changes

nrow(seas[nuc.div.is.seasonal == TRUE & abund.is.seasonal ==TRUE]) / nrow(seas) # 24% of the abundant-enough genomes had both seasonal diversity and abundance patterns

nrow(seas[blooms.are == "Less diverse"]) / nrow(seas[nuc.div.is.seasonal == TRUE & abund.is.seasonal ==TRUE]) # 0.2054795 % had more clonal blooms

nrow(seas[blooms.are == "More diverse"]) / nrow(seas[nuc.div.is.seasonal == TRUE & abund.is.seasonal ==TRUE]) # 0.1945205 % had more diverse blooms

nrow(snv[Classified.LT.Change == "gradual" | Classified.LT.Change == "step" | Classified.LT.Change == "disturbance"])/nrow(snv) * 100 # 20.91255 % had long-term patterns

nrow(snv[Classified.LT.Change == "step" | Classified.LT.Change == "disturbance"]) / nrow(snv[Classified.LT.Change == "gradual" | Classified.LT.Change == "step" | Classified.LT.Change == "disturbance"]) * 100 # 65% of LT petterns were abrupt change

nrow(snv[Classified.LT.Change == "gradual"] )/nrow(snv) * 100 # 7.224335 % had gradual change

nrow(snv[Classified.LT.Change == "step" ])/nrow(snv) * 100 # 6.08365 had step changes

nrow(snv[Classified.LT.Change == "disturbance"])/nrow(snv) * 100 # 7.604563 had step changes

nrow(tax[family == "Nanopelagicaceae"]) # 127 acI genomes recovered

sum(tax[family == "Nanopelagicaceae", sum.abund]) / sum(tax[ ,sum.abund]) * 100 # acI accounts for 19.00661 % of total mapped reads

# ---- pull out barplot table ----

seas <- merge(x = seas, y = tax, by = "genome", all = F)

get.plotting.matrix <- function(xaxis.choice = "phylum",
                                yaxis.choice = "num.genomes",
                                num.bars = 6,
                                as.percent = TRUE){
  ss <- seas[ , .(num.genomes = .N, mean.abund = mean(mean.abund)), by = c(xaxis.choice, "abund.is.seasonal", "nuc.div.is.seasonal")]
  form <- paste("abund.is.seasonal + nuc.div.is.seasonal ~", xaxis.choice)
  ss <- dcast(data = ss, formula = form, value.var = "num.genomes", fill = 0)
  ss[abund.is.seasonal == FALSE & nuc.div.is.seasonal == TRUE , label := "Div."]
  ss[abund.is.seasonal == TRUE & nuc.div.is.seasonal == TRUE , label := "Both"]
  ss[abund.is.seasonal == TRUE & nuc.div.is.seasonal == FALSE , label := "Abund"]
  ss[abund.is.seasonal == FALSE & nuc.div.is.seasonal == FALSE , label := "Neither"]
  ss <- as.matrix(x = ss[ ,-c(1:2)], rownames = "label")
  ss <- ss[c(1, 3,4,2), ]
  col.order <- order(colSums(ss), decreasing = T)
  ss <- ss[ ,col.order]
  others <- ss[ ,(num.bars + 1):ncol(ss)]
  ss <- cbind(ss[ ,1:num.bars], "Other" = rowSums(others))
  ss <- as.matrix(ss)
  if (as.percent){
     tots <- colSums(ss)
     ss <- t(ss)
     ss <- ss / tots * 100
     ss <- t(ss)
     # colSums(ss) # check 100
  }
  return(ss)
}

ss <- get.plotting.matrix(xaxis.choice = "phylum", yaxis.choice = "num.genomes", num.bars = 6, as.percent = TRUE)

# ---- overall how seasonal ----

par(mar = c(4,3.5,.5,7.5))
y.lim <- max(colSums(ss))
bar.spots <- barplot(height = ss, legend.text = FALSE, col = c("grey80","#e79493","#9865c0","#a9adff"), axes = F, ann = F, names.arg = rep("",ncol(ss)))
text(x = bar.spots, y = -(y.lim / 50), labels = colnames(ss), xpd = NA, srt = 30, adj = 1)
axis(side = 2, at = c(0,y.lim), labels = F, lwd = 1, lwd.ticks = 0, line = -.5)
axis(side = 2, labels = F, line = -.5, lwd = 0, lwd.ticks = 1, tck = -.02)
axis(side = 2, labels = T, line = -.75, lwd = 0, lwd.ticks = 0, las = 2)
mtext(text = "Percent Genomes", side = 2, line = 2.5, outer = F)
y.legend.locs <- (y.lim * .8) - (y.lim / c(15))*(1:5)
rect(xleft = max(bar.spots) + 1, xright = max(bar.spots) + 1.5, ytop = y.legend.locs[-1] - (y.lim / 50), ybottom = y.legend.locs[-1] + (y.lim / 50), xpd = NA, col = c("#a9adff","#9865c0","#e79493","grey80"))
text(x = max(bar.spots) + c(1,rep(1.75,4)), y = y.legend.locs, labels = c("Seasonality in:","Nucleotide\ndiversity","Both","Abundance","Neither"), xpd = NA, adj = 0)

# ---- overall how bloomswq







# get colors: 
plot(x = 1:3,y = 1:3, type = "n")
points(x = 1:2, y = 1:2, pch = 19, cex = 5, col = adjustcolor("red3",.5))
points(x = 2:3, y = 2:3, pch = 19, cex = 5, col = adjustcolor("blue",.4))






