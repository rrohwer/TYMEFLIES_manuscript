# RRR

# seriously do I really need to go back and remove short contigs?
# look at all the bins from sample 1, are short contigs even there??

library(phylotools)
library(tidyr)

options(scipen = 9999)

plot.folder <- "figures/2022-11-16_short_contig_contributions"

all.bins <- readRDS(file = "data/2022-11-15_bin_stats/all_bins.rds")
all.bins <- all.bins[all.bins$tymeflies.name == "ME2000-03-15pf_3300044539", ]

bin.files <- list.files("data/2022-08-02_example_bins/ME2000-03-15pf_3300044539", full.names = T)
bin.names <- list.files("data/2022-08-02_example_bins/ME2000-03-15pf_3300044539", full.names = F)
bin.names <- sub("\\.fa","",bin.names)
bin.names <- sub("\\.","",bin.names)

my.bins <- list(NULL)
for (b in 1:length(bin.files)){
  my.bins[[b]] <- read.fasta(file = bin.files[b])
  my.bins[[b]]$bin <- bin.names[b]
}

bin.table <- do.call(what = rbind, args = my.bins)
head(bin.table, n = 2)

bin.table$length <- nchar(bin.table$seq.text)

summary(bin.table$length)
bin.table$length

hist(bin.table$length)

raw.bin.table <- bin.table
raw.all.bins <- all.bins

# ---- 4000 bp cutoff ----

bin.table <- raw.bin.table
all.bins <- raw.all.bins

bin.table$is.small <- raw.bin.table$length < 4000

long.bp <- sum(bin.table$length[!bin.table$is.small])
short.bp <- sum(bin.table$length[bin.table$is.small])

pie(x = c(long.bp,short.bp), labels = c("long","short"))

colnames(bin.table)
x <- pivot_wider(data = bin.table, id_cols = "bin", names_from = "is.small", values_from = "length", values_fn = sum)
x <- as.data.frame(x)
colnames(x) <- c("bin","long","short")

barplot(height = t(as.matrix(x[ ,-1])), beside = F, names.arg = x[ ,1], legend = T, las = 2, cex.names = .5)

head(x)
colnames(all.bins)
all.bins <- merge(x = x, y = all.bins, by.x = "bin", by.y = "bin.id")
all.bins$color <- "grey"
all.bins$color[all.bins$completeness >= 50 & all.bins$contamination < 10] <- "orange"
all.bins$color[all.bins$completeness >= 90 & all.bins$contamination < 5] <- "red"

low.contam <- all.bins$contamination < 1
HQ <- all.bins$completeness >= 90 & all.bins$contamination < 5
MQ <- all.bins$completeness >= 50 & all.bins$contamination < 10 & !HQ
acI <- all.bins$family == "f__Nanopelagicaceae"
acI[is.na(acI)] <- FALSE
was.winner <- all.bins$winner == TRUE

png(filename = file.path(plot.folder,"sample_1_all_bins.png"), width = 20, height = 4, units = "in", res = 300)
par(mar = c(6,4,1,0))
barplot(height = t(as.matrix(all.bins[ ,c("long","short")])), beside = F, names.arg = all.bins$family, legend = T, las = 2, cex.names = .5, axes = F)
axis(side = 2, las = 2, line = -3)
mtext("all sample_1 bins", cex = 1.5, line = -1)
mtext("Contribution to Bin (bp)", side = 2, line = 2.5, cex = 1.5)
dev.off()

png(filename = file.path(plot.folder,"sample_1_MQ_bins.png"), width = 10, height = 4, units = "in", res = 300)
par(mar = c(7.5,5,1,0))
barplot(height = t(as.matrix(all.bins[MQ,c("long","short")])), beside = F, names.arg = all.bins$family[MQ], legend = T, las = 2, cex.names = .7, axes = F)
mtext("sample_1 MQ bins")
axis(side = 2, las = 2, line = -1)
mtext("Contribution to Bin (bp)", side = 2, line = 4, cex = 1)
dev.off()

png(filename = file.path(plot.folder,"sample_1_HQ_bins.png"), width = 6.5, height = 4, units = "in", res = 300)
par(mar = c(9,5,1,0))
barplot(height = t(as.matrix(all.bins[HQ,c("long","short")])), beside = F, names.arg = all.bins$family[HQ], legend = T, las = 2, cex.names = .9, axes = F)
mtext("sample_1 HQ bins")
axis(side = 2, las = 2, line = -1, cex = .9)
mtext("Contribution to Bin (bp)", side = 2, line = 4, cex = 1)
dev.off()

png(filename = file.path(plot.folder,"sample_1_contam_less_1_bins.png"), width = 15, height = 4, units = "in", res = 300)
par(mar = c(6,4,1,0))
barplot(height = t(as.matrix(all.bins[low.contam,c("long","short")])), beside = F, names.arg = all.bins$family[low.contam], legend = T, las = 2, cex.names = .5, axes = F)
mtext("sample_1 contam < 1% bins")
axis(side = 2, las = 2, line = -2, cex = .9)
mtext("Contribution to Bin (bp)", side = 2, line = 3, cex = 1)
dev.off()

png(filename = file.path(plot.folder,"sample_1_acI_bins.png"), width = 6.5, height = 4, units = "in", res = 300)
par(mar = c(5.5,5,1,0))
barplot(height = t(as.matrix(all.bins[acI,c("long","short")])), beside = F, names.arg = all.bins$family[acI], legend = T, las = 2, cex.names = .5, axes = F)
mtext("sample_1 acI bins")
axis(side = 2, las = 2, line = -.5, cex = .9)
mtext("Contribution to Bin (bp)", side = 2, line = 4, cex = 1)
dev.off()

png(filename = file.path(plot.folder,"sample_1_drep_winner_bins.png"), width = 6.5, height = 4, units = "in", res = 300)
par(mar = c(5.5,5,2,0))
barplot(height = t(as.matrix(all.bins[was.winner,c("long","short")]))[ ,c(1,2,3,5,4)], beside = F, names.arg = all.bins$family[was.winner][c(1,2,3,5,4)], legend = T, las = 2, cex.names = .5, axes = F)
mtext("sample_1 drep winner bins", line = 1)
axis(side = 2, las = 2, line = -.5, cex = .9)
mtext("Contribution to Bin (bp)", side = 2, line = 4, cex = 1)
dev.off()

# ---- 3000 bp cutoff ----

bin.table <- raw.bin.table
all.bins <- raw.all.bins

bin.table$is.small <- raw.bin.table$length < 3000

long.bp <- sum(bin.table$length[!bin.table$is.small])
short.bp <- sum(bin.table$length[bin.table$is.small])

pie(x = c(long.bp,short.bp), labels = c("long","short"))

colnames(bin.table)
x <- pivot_wider(data = bin.table, id_cols = "bin", names_from = "is.small", values_from = "length", values_fn = sum)
x <- as.data.frame(x)
colnames(x) <- c("bin","long","short")

barplot(height = t(as.matrix(x[ ,-1])), beside = F, names.arg = x[ ,1], legend = T, las = 2, cex.names = .5)

head(x)
colnames(all.bins)
all.bins <- merge(x = x, y = all.bins, by.x = "bin", by.y = "bin.id")
all.bins$color <- "grey"
all.bins$color[all.bins$completeness >= 50 & all.bins$contamination < 10] <- "orange"
all.bins$color[all.bins$completeness >= 90 & all.bins$contamination < 5] <- "red"
      
low.contam <- all.bins$contamination < 1
HQ <- all.bins$completeness >= 90 & all.bins$contamination < 5
MQ <- all.bins$completeness >= 50 & all.bins$contamination < 10 & !HQ
acI <- all.bins$family == "f__Nanopelagicaceae"
acI[is.na(acI)] <- FALSE
was.winner <- all.bins$winner == TRUE
    
png(filename = file.path(plot.folder,"sample_1_all_bins-3000.png"), width = 20, height = 4, units = "in", res = 300)
par(mar = c(6,4,1,0))
barplot(height = t(as.matrix(all.bins[ ,c("long","short")])), beside = F, names.arg = all.bins$family, legend = T, las = 2, cex.names = .5, axes = F)
axis(side = 2, las = 2, line = -3)
mtext("all sample_1 bins", cex = 1.5, line = -1)
mtext("Contribution to Bin (bp)", side = 2, line = 2.5, cex = 1.5)
dev.off()

png(filename = file.path(plot.folder,"sample_1_MQ_bins-3000.png"), width = 10, height = 4, units = "in", res = 300)
par(mar = c(7.5,5,1,0))
barplot(height = t(as.matrix(all.bins[MQ,c("long","short")])), beside = F, names.arg = all.bins$family[MQ], legend = T, las = 2, cex.names = .7, axes = F)
mtext("sample_1 MQ bins")
axis(side = 2, las = 2, line = -1)
mtext("Contribution to Bin (bp)", side = 2, line = 4, cex = 1)
dev.off()

png(filename = file.path(plot.folder,"sample_1_HQ_bins-3000.png"), width = 6.5, height = 4, units = "in", res = 300)
par(mar = c(9,5,1,0))
barplot(height = t(as.matrix(all.bins[HQ,c("long","short")])), beside = F, names.arg = all.bins$family[HQ], legend = T, las = 2, cex.names = .9, axes = F)
mtext("sample_1 HQ bins")
axis(side = 2, las = 2, line = -1, cex = .9)
mtext("Contribution to Bin (bp)", side = 2, line = 4, cex = 1)
dev.off()

png(filename = file.path(plot.folder,"sample_1_contam_less_1_bins-3000.png"), width = 15, height = 4, units = "in", res = 300)
par(mar = c(6,4,1,0))
barplot(height = t(as.matrix(all.bins[low.contam,c("long","short")])), beside = F, names.arg = all.bins$family[low.contam], legend = T, las = 2, cex.names = .5, axes = F)
mtext("sample_1 contam < 1% bins")
axis(side = 2, las = 2, line = -2, cex = .9)
mtext("Contribution to Bin (bp)", side = 2, line = 3, cex = 1)
dev.off()

png(filename = file.path(plot.folder,"sample_1_acI_bins-3000.png"), width = 6.5, height = 4, units = "in", res = 300)
par(mar = c(5.5,5,1,0))
barplot(height = t(as.matrix(all.bins[acI,c("long","short")])), beside = F, names.arg = all.bins$family[acI], legend = T, las = 2, cex.names = .5, axes = F)
mtext("sample_1 acI bins")
axis(side = 2, las = 2, line = -.5, cex = .9)
mtext("Contribution to Bin (bp)", side = 2, line = 4, cex = 1)
dev.off()

png(filename = file.path(plot.folder,"sample_1_drep_winner_bins-3000.png"), width = 6.5, height = 4, units = "in", res = 300)
par(mar = c(5.5,5,2,0))
barplot(height = t(as.matrix(all.bins[was.winner,c("long","short")]))[ ,c(1,2,3,5,4)], beside = F, names.arg = all.bins$family[was.winner][c(1,2,3,5,4)], legend = T, las = 2, cex.names = .5, axes = F)
mtext("sample_1 drep winner bins", line = 1)
axis(side = 2, las = 2, line = -.5, cex = .9)
mtext("Contribution to Bin (bp)", side = 2, line = 4, cex = 1)
dev.off()

# I'm concerned about throwing so much of the bin contents out. 

