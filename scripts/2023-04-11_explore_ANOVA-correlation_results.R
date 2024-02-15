# RRR

# I re-did the correlations on only the nuc div with temporal patterns
# and with more lenient inclusion cutoffs otherwise
# now re-look at the new results
# use baseR barplots because ggplot barplots suuuuck

options(scipen = 9999)

library(data.table)
library(lubridate)

s.cor <- fread(file = "data/2023-04-09_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-SEASONAL.tsv")
a.cor <- fread(file = "data/2023-04-09_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-ANNUAL.tsv")

tax <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
tax <- as.data.table(tax)

rel.abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
colnames(rel.abund)[3] <- "abund.perc"

cor.cutoff <- .3 # debate later!!

num.bars <- 20

# ---- functions ----

summarize.by.tax.level <- function(my.cor, tax.level, num.bars = NA){
  my.cor <- my.cor[ , .(count = .N, abund = sum(sum.abund)), by = c(tax.level, "match")]
  if (!is.na(num.bars)){
    if (nrow(unique(my.cor[ ,..tax.level])) > num.bars){
      my.cor <- my.cor[order(-abund)]
      my.genomes <- unique(my.cor[ ,..tax.level])
      keep <- my.genomes[1:num.bars]
      other <- my.genomes[(num.bars + 1):nrow(my.genomes)]
      keep <- merge(x = my.cor, y = keep, by = tax.level, all.x = F, all.y = T)
      other <- merge(x = my.cor, y = other, by = tax.level, all.x = F, all.y = T)
      other <- other[ ,1 := "Other"]
      other <- other[ , .(count = sum(count), abund = sum(abund)), by = c(tax.level,"match")]
      my.cor <- rbind(keep, other)
    }
  }
  num.genomes <- my.cor[ , .(num.genomes = sum(count)), by = tax.level]
  num.genomes <- num.genomes[order(-num.genomes)]
  my.cor <- merge(x = my.cor, y = num.genomes, by = tax.level, all = T)
  return(my.cor)
}

get.percent.summary.table <- function(my.tax, my.column, tax.level){ 
  # my.column = "abund" or "count"
  x <- dcast(data = my.tax, formula = paste(tax.level, "~ match"), value.var = my.column)
  x <- as.matrix(x, rownames = T)
  x[is.na(x)] <- 0
  x <- x[order(rowSums(x), decreasing = T), ]
  x <- x / rowSums(x) * 100
  x <- t(x)
  x <- x[c(3,1,2), ]
  return(x)
}

get.absolute.summary.table <- function(my.tax, my.column, tax.level){ 
  # my.column = "abund" or "count"
  x <- dcast(data = my.tax, formula = paste(tax.level, "~ match"), value.var = my.column)
  x <- as.matrix(x, rownames = T)
  x[is.na(x)] <- 0
  x <- x[order(rowSums(x), decreasing = T), ]
  x <- t(x)
  x <- x[c(3,1,2), ]
  return(x)
}

make.barplot <- function(my.mat, my.tax, include.other = FALSE){
  # need matrix for the barplot and order of bars
  # need tax table for number genomes information
  if(!include.other){
    index <- which(colnames(my.mat) == "Other")
    if(length(index) > 0){
      my.mat <- my.mat[ ,-index]
    }
  }
  bar.spots <- barplot(height = my.mat, col = c("grey70", "gold2", "aquamarine3"), names.arg = rep("",ncol(my.mat)), axes = F)
  plot.info <- data.table("loc" = bar.spots, "tax" = colnames(my.mat), "max" = colSums(my.mat))
  my.tax <- my.tax[ ,c(..t,"num.genomes")]
  my.tax <- merge(x = plot.info, y = my.tax, by.x = "tax", by.y = t, all.x = TRUE, all.y = FALSE )
  my.tax <- unique(my.tax)
  text(x = my.tax$loc, y = -max(colSums(my.mat)) / 20, labels = my.tax$tax, srt = 90, xpd = NA, adj = 1, cex = .7)
  text(x = my.tax$loc, y = my.tax$max + (max(colSums(my.mat)) / 20), labels = my.tax$num.genomes, xpd = NA, adj = .5, cex = .7)
  axis(side = 2, cex.axis = .7, las = 2, line = -.25, labels = F, tck = -.025)
  axis(side = 2, cex.axis = .7, las = 2, line = -.75, labels = T, lwd = 0)
}

# ---- clean up taxonomy names ----

rel.abund <- rel.abund[ ,.(sum.abund = sum(abund.perc)), by = Genome]
tax <- merge(x = rel.abund, y = tax, by.x = "Genome", by.y = "bin.full.name", all = F)

tax[ ,`:=`(phylum = sub(pattern = "p__", replacement = "", x = phylum), 
           class = sub(pattern = "c__", replacement = "", x = class),
           order = sub(pattern = "o__", replacement = "", x = order),
           family = sub(pattern = "f__", replacement = "", x = family),
           genus = sub(pattern = "g__", replacement = "", x = genus),
           species = sub(pattern = "s__", replacement = "", x = species))]

tax[is.na(phylum) | phylum == "", phylum := "Unclassified"]
tax[is.na(class) | class == "", class := paste0("U.",phylum)]
tax[is.na(order) | order == "", order := paste0("U.",class)]
tax[is.na(family) | family == "", family := paste0("U.",order)]
tax[is.na(genus) | genus == "", genus := paste0("U.",family)]
tax[is.na(species) | species == "", species := paste0("U.",genus)]

tax[ , `:=`(class = sub(pattern = "^(U\\.)+", replacement = "U\\.", x = class),
            order = sub(pattern = "^(U\\.)+", replacement = "U\\.", x = order),
            family = sub(pattern = "^(U\\.)+", replacement = "U\\.", x = family),
            genus = sub(pattern = "^(U\\.)+", replacement = "U\\.", x = genus),
            species = sub(pattern = "^(U\\.)+", replacement = "U\\.", x = species))]

# ---- parse correlation data ----

s.cor <- merge(x = s.cor, y = tax, by.x = "genome", by.y = "Genome", all.x = TRUE, all.y = FALSE)
a.cor <- merge(x = a.cor, y = tax, by.x = "genome", by.y = "Genome", all.x = TRUE, all.y = FALSE)

s.cor[cor >= cor.cutoff, match := "match"]
s.cor[cor <= -cor.cutoff, match := "offset"]
s.cor[cor > -cor.cutoff & cor < cor.cutoff, match := "other"]
s.cor[ , match := factor(match, levels = c("match","offset","other"))]

a.cor[cor >= cor.cutoff, match := "match"]
a.cor[cor <= -cor.cutoff, match := "offset"]
a.cor[cor > -cor.cutoff & cor < cor.cutoff, match := "other"]
a.cor[ , match := factor(match, levels = c("match","offset","other"))]

# ---- summarize by taxa level and plot ----

tax.levels <- c("phylum","class","order","family","genus","species")
for (t in tax.levels){
  
  a.tax <- summarize.by.tax.level(my.cor = a.cor, tax.level = t, num.bars = num.bars)
  s.tax <- summarize.by.tax.level(my.cor = s.cor, tax.level = t, num.bars = num.bars)
  
  a.abund.perc <- get.percent.summary.table(my.tax = a.tax, my.column = "abund", tax.level = t) 
  s.abund.perc <- get.percent.summary.table(my.tax = s.tax, my.column = "abund", tax.level = t) 
  
  a.count.perc <- get.percent.summary.table(my.tax = a.tax, my.column = "count", tax.level = t) 
  s.count.perc <- get.percent.summary.table(my.tax = s.tax, my.column = "count", tax.level = t) 
  
  a.abund.abs <- get.absolute.summary.table(my.tax = a.tax, my.column = "abund", tax.level = t)
  s.abund.abs <- get.absolute.summary.table(my.tax = s.tax, my.column = "abund", tax.level = t)
  
  a.count.abs <- get.absolute.summary.table(my.tax = a.tax, my.column = "count", tax.level = t)
  s.count.abs <- get.absolute.summary.table(my.tax = s.tax, my.column = "count", tax.level = t)
  # ----
 
  pdf(file = paste0("figures/2023-04-11_div-abund_correlation_barplots/by-abundance_percent-of-taxon/",which(tax.levels == t),"-",t,".pdf"), width = 6.5, height = 4)
  par(mar = c(7,1,2.5,0), oma = c(.1,2.5,.1,0), mfrow = c(1,2))
  make.barplot(my.mat = a.abund.perc, my.tax = a.tax)
  mtext(text = paste0("Percent of ",t,"\n(scaled by sum abundance)"), side = 2, line = 1.5, cex = .9)
  mtext(text = c("Long-Term"), side = 3, outer = T, at = c(.3), col = c("black"), cex = 1, line = -1)
  make.barplot(my.mat = s.abund.perc, my.tax = s.tax)
  mtext(text = c("Seasonal"), side = 3, outer = T, at = c(.8), col = c("black"), cex = 1, line = -1)
  mtext(text = c("match","offset","other"), side = 3, outer = T, at = c(-.06), col = c("aquamarine3", "gold2","grey70"), cex = 1, line = c(-16,-17,-18), adj = 0)
  dev.off()
  
  pdf(file = paste0("figures/2023-04-11_div-abund_correlation_barplots/by-abundance_absolute_abund/",which(tax.levels == t),"-",t,".pdf"), width = 6.5, height = 4)
  par(mar = c(7,1,2.5,0), oma = c(.1,2.5,.1,0), mfrow = c(1,2))
  make.barplot(my.mat = a.abund.abs, my.tax = a.tax)
  mtext(text = paste0("Abundance of ",t,"\n(sum percent abundance)"), side = 2, line = 1.5, cex = .9)
  mtext(text = c("Long-Term"), side = 3, outer = T, at = c(.3), col = c("black"), cex = 1, line = -1)
  make.barplot(my.mat = s.abund.abs, my.tax = s.tax)
  mtext(text = c("Seasonal"), side = 3, outer = T, at = c(.8), col = c("black"), cex = 1, line = -1)
  mtext(text = c("match","offset","other"), side = 3, outer = T, at = c(-.06), col = c("aquamarine3", "gold2","grey70"), cex = 1, line = c(-16,-17,-18), adj = 0)
  dev.off()
  
  pdf(file = paste0("figures/2023-04-11_div-abund_correlation_barplots/by-count_percent-of-taxon/",which(tax.levels == t),"-",t,".pdf"), width = 6.5, height = 4)
  par(mar = c(7,1,2.5,0), oma = c(.1,2.5,.1,0), mfrow = c(1,2))
  make.barplot(my.mat = a.count.perc, my.tax = a.tax)
  mtext(text = paste0("Percent of ",t,"\n(scaled by sum total)"), side = 2, line = 1.5, cex = .9)
  mtext(text = c("Long-Term"), side = 3, outer = T, at = c(.3), col = c("black"), cex = 1, line = -1)
  make.barplot(my.mat = s.count.perc, my.tax = s.tax)
  mtext(text = c("Seasonal"), side = 3, outer = T, at = c(.8), col = c("black"), cex = 1, line = -1)
  mtext(text = c("match","offset","other"), side = 3, outer = T, at = c(-.06), col = c("aquamarine3", "gold2","grey70"), cex = 1, line = c(-16,-17,-18), adj = 0)
  dev.off()
  
  pdf(file = paste0("figures/2023-04-11_div-abund_correlation_barplots/by-count_absolute_count/",which(tax.levels == t),"-",t,".pdf"), width = 6.5, height = 4)
  par(mar = c(7,1,2.5,0), oma = c(.1,2.5,.1,0), mfrow = c(1,2))
  make.barplot(my.mat = a.count.abs, my.tax = a.tax)
  mtext(text = paste0("Count of ",t,"\n(total count)"), side = 2, line = 1.5, cex = .9)
  mtext(text = c("Long-Term"), side = 3, outer = T, at = c(.3), col = c("black"), cex = 1, line = -1)
  make.barplot(my.mat = s.count.abs, my.tax = s.tax)
  mtext(text = c("Seasonal"), side = 3, outer = T, at = c(.8), col = c("black"), cex = 1, line = -1)
  mtext(text = c("match","offset","other"), side = 3, outer = T, at = c(-.06), col = c("aquamarine3", "gold2","grey70"), cex = 1, line = c(-16,-17,-18), adj = 0)
  dev.off()
  
  pdf(file = paste0("figures/2023-04-11_div-abund_correlation_barplots/by-abundance_percent-of-taxon/",which(tax.levels == t),"-",t,"_with_other.pdf"), width = 6.5, height = 4)
  par(mar = c(7,1,2.5,0), oma = c(.1,2.5,.1,0), mfrow = c(1,2))
  make.barplot(my.mat = a.abund.perc, my.tax = a.tax, include.other = TRUE)
  mtext(text = paste0("Percent of ",t,"\n(scaled by sum abundance)"), side = 2, line = 1.5, cex = .9)
  mtext(text = c("Long-Term"), side = 3, outer = T, at = c(.3), col = c("black"), cex = 1, line = -1)
  make.barplot(my.mat = s.abund.perc, my.tax = s.tax, include.other = TRUE)
  mtext(text = c("Seasonal"), side = 3, outer = T, at = c(.8), col = c("black"), cex = 1, line = -1)
  mtext(text = c("match","offset","other"), side = 3, outer = T, at = c(-.06), col = c("aquamarine3", "gold2","grey70"), cex = 1, line = c(-16,-17,-18), adj = 0)
  dev.off()
  
  pdf(file = paste0("figures/2023-04-11_div-abund_correlation_barplots/by-abundance_absolute_abund/",which(tax.levels == t),"-",t,"_with_other.pdf"), width = 6.5, height = 4)
  par(mar = c(7,1,2.5,0), oma = c(.1,2.5,.1,0), mfrow = c(1,2))
  make.barplot(my.mat = a.abund.abs, my.tax = a.tax, include.other = TRUE)
  mtext(text = paste0("Abundance of ",t,"\n(sum percent abundance)"), side = 2, line = 1.5, cex = .9)
  mtext(text = c("Long-Term"), side = 3, outer = T, at = c(.3), col = c("black"), cex = 1, line = -1)
  make.barplot(my.mat = s.abund.abs, my.tax = s.tax, include.other = TRUE)
  mtext(text = c("Seasonal"), side = 3, outer = T, at = c(.8), col = c("black"), cex = 1, line = -1)
  mtext(text = c("match","offset","other"), side = 3, outer = T, at = c(-.06), col = c("aquamarine3", "gold2","grey70"), cex = 1, line = c(-16,-17,-18), adj = 0)
  dev.off()
  
  pdf(file = paste0("figures/2023-04-11_div-abund_correlation_barplots/by-count_percent-of-taxon/",which(tax.levels == t),"-",t,"_with_other.pdf"), width = 6.5, height = 4)
  par(mar = c(7,1,2.5,0), oma = c(.1,2.5,.1,0), mfrow = c(1,2))
  make.barplot(my.mat = a.count.perc, my.tax = a.tax, include.other = TRUE)
  mtext(text = paste0("Percent of ",t,"\n(scaled by sum total)"), side = 2, line = 1.5, cex = .9)
  mtext(text = c("Long-Term"), side = 3, outer = T, at = c(.3), col = c("black"), cex = 1, line = -1)
  make.barplot(my.mat = s.count.perc, my.tax = s.tax, include.other = TRUE)
  mtext(text = c("Seasonal"), side = 3, outer = T, at = c(.8), col = c("black"), cex = 1, line = -1)
  mtext(text = c("match","offset","other"), side = 3, outer = T, at = c(-.06), col = c("aquamarine3", "gold2","grey70"), cex = 1, line = c(-16,-17,-18), adj = 0)
  dev.off()
  
  pdf(file = paste0("figures/2023-04-11_div-abund_correlation_barplots/by-count_absolute_count/",which(tax.levels == t),"-",t,"_with_other.pdf"), width = 6.5, height = 4)
  par(mar = c(7,1,2.5,0), oma = c(.1,2.5,.1,0), mfrow = c(1,2))
  make.barplot(my.mat = a.count.abs, my.tax = a.tax, include.other = TRUE)
  mtext(text = paste0("Count of ",t,"\n(total count)"), side = 2, line = 1.5, cex = .9)
  mtext(text = c("Long-Term"), side = 3, outer = T, at = c(.3), col = c("black"), cex = 1, line = -1)
  make.barplot(my.mat = s.count.abs, my.tax = s.tax, include.other = TRUE)
  mtext(text = c("Seasonal"), side = 3, outer = T, at = c(.8), col = c("black"), cex = 1, line = -1)
  mtext(text = c("match","offset","other"), side = 3, outer = T, at = c(-.06), col = c("aquamarine3", "gold2","grey70"), cex = 1, line = c(-16,-17,-18), adj = 0)
  dev.off()
  
  # ----
}












