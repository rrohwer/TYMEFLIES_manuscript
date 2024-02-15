# RRR
# make simplified genes under selection plots for ESA
# Don't include the COG categories, just do totals

# This time I normalized by number of samples- this did dampen some of the "trends"
# But now I'm wondering, in the future do I also need to normalize by abundance? b/c higher detection limit when more abundant? 
# Or maybe just exclude samples below X coverage. ooh yes this.

library(data.table)
library(lubridate)

acI.B <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz")
genes <- acI.B
my.genome <- "ME2011-09-21_3300043464_group3_bin69"

acI.C <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2016-07-20_3300033996_group7_bin32_gene_info_combined_MK.tsv.gz")
genes <- acI.C
my.genome <- "ME2016-07-20_3300033996_group7_bin32"

acI.A <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-04_3300044729_group3_bin142_gene_info_combined_MK.tsv.gz")
genes <- acI.A
my.genome <- "ME2011-09-04_3300044729_group3_bin142"
# ---- parse data ----

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,sample := sub(pattern = "\\.IS_gene_info\\.tsv", replacement = "", x = sample)]
genes[ ,sample := sub(pattern = ".gz", replacement = "", x = sample)]

pos <- genes[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman), .(.N), by = .(sample, date, season, year)]
neg <- genes[MK.p.val >= .05 & mcdonald.kreitman > 1 & is.finite(mcdonald.kreitman), .(.N), by = .(sample, date, season, year)] # AHHH BUG BUG BUG, p-value should be LESS than 0.05!!!!
pos <- pos[ ,.(N = mean(N)), by = .(date, season, year)]
neg <- neg[ ,.(N = mean(N)), by = .(date, season, year)]

sel.summary <- merge(x = pos, y = neg, by = c("date", "season", "year"), all = T, suffixes = c("pos","neg"))
sel.summary[is.na(Npos), Npos := 0]
sel.summary[is.na(Nneg), Nneg := 0]

year.sum <- sel.summary[ ,.(.N, Npos = sum(Npos), Nneg = sum(Nneg)), by = year]
year.sum[ ,`:=`(Npos.per.sample = Npos / N,
                Nneg.per.sample = Nneg / N)]

season.sum <- sel.summary[ ,.(.N, Npos = sum(Npos), Nneg = sum(Nneg)), by = season]
season.sum[ ,`:=`(Npos.per.sample = Npos / N,
                Nneg.per.sample = Nneg / N)]

season.order <- data.table("order" = 1:6, 
                           season = c("Spring","Clearwater","Early.Summer","Late.Summer","Fall","Ice.On"),
                           label = c("Spring","Clearwater","Early\nSummer","Late\nSummer","Fall","Ice-On"))
season.sum <- merge(x = season.sum, y = season.order, by = "season")
season.sum <- season.sum[order(order)]

# ---- make plots ----

# for acI-A ----
pdf(file = "figures/2023-08-05_ESA_plots/gene_selection_acI-A_seasons.pdf", width = 6.5, height = 3)
bar.spots <- barplot(height = season.sum$Npos.per.sample, beside = T, plot = F, space = .75)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(season.sum$Nneg.per.sample, na.rm = T), max(season.sum$Npos.per.sample, na.rm = T))
par(mar = c(2.5,6,2,1))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = season.sum$Npos.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = "#b89996", add = T, border = NA, space = .75)
barplot(height = -season.sum$Nneg.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = "#7592b2", add = T, border = NA, space = .75)
segments(x0 = .3, x1 = 11, y0 = 0, y1 = 0, xpd = T)
text(x = bar.spots, y = -28, labels = season.sum$label, cex = 1, xpd = NA)
y.lim
  ticks <- axis(side = 2, line = 1, labels = F, tck = -.02, at = seq.int(-25,5,5), xpd = T)
axis(side = 2, line = .75, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = 1, xpd = T)
mtext(text = "Genes Under Selection (ave / sample)", side = 2, line = 5, cex = 1)
mtext(text = "Positive", side = 2, line = 3.5, at = 1, cex = 1, adj = 0)
mtext(text = "Purifying", side = 2, line = 3.5, at = -1, cex = 1, adj = 1)
x.lim
arrows(x0 = -.3, y0 = 1, x1 = -.3, y1 = 9.5, xpd = T, length = .075)
arrows(x0 = -.3, y0 = -1, x1 = -.3, y1 = -11, xpd = T, length = .075)
dev.off()
# ----

# for acI-C ----
pdf(file = "figures/2023-08-05_ESA_plots/gene_selection_acI-C_seasons.pdf", width = 6.5, height = 3)
bar.spots <- barplot(height = season.sum$Npos.per.sample, beside = T, plot = F, space = .75)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(season.sum$Nneg.per.sample, na.rm = T), max(season.sum$Npos.per.sample, na.rm = T))
par(mar = c(2.5,6,2,1))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = season.sum$Npos.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = "#b89996", add = T, border = NA, space = .75)
barplot(height = -season.sum$Nneg.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = "#7592b2", add = T, border = NA, space = .75)
segments(x0 = .3, x1 = 11, y0 = 0, y1 = 0, xpd = T)
text(x = bar.spots, y = -45, labels = season.sum$label, cex = 1, xpd = NA)
y.lim
ticks <- axis(side = 2, line = 1, labels = F, tck = -.02, at = seq.int(-40,10,10), xpd = T)
axis(side = 2, line = .75, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = 1, xpd = T)
mtext(text = "Genes Under Selection (ave / sample)", side = 2, line = 5, cex = 1)
mtext(text = "Positive", side = 2, line = 3.5, at = 2, cex = 1, adj = 0)
mtext(text = "Purifying", side = 2, line = 3.5, at = -2, cex = 1, adj = 1)
x.lim
arrows(x0 = -.3, y0 = 2, x1 = -.3, y1 = 17, xpd = T, length = .075)
arrows(x0 = -.3, y0 = -2, x1 = -.3, y1 = -19, xpd = T, length = .075)
dev.off()
# ----

# for acI-B ----
pdf(file = "figures/2023-08-05_ESA_plots/gene_selection_acI-B_long-term.pdf", width = 6.5, height = 3)
bar.spots <- barplot(height = year.sum$Npos.per.sample, beside = T, plot = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(year.sum$Nneg.per.sample, na.rm = T), max(year.sum$Npos.per.sample, na.rm = T))
par(mar = c(2.5,5,2,0))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = year.sum$Npos.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = "#b89996", add = T, border = NA)
barplot(height = -year.sum$Nneg.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = "#7592b2", add = T, border = NA)
abline(h = 0)
text(x = bar.spots[seq.int(1,20,2)], y = -52, labels = year.sum$year[seq.int(1,20,2)], cex = 1, xpd = NA, srt = 90)
y.lim
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02, at = seq.int(-50,10,10), xpd = T)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = 1, xpd = T)
mtext(text = "Genes Under Selection (ave / sample)", side = 2, line = 4, cex = 1, at = -20)
mtext(text = "Positive", side = 2, line = 2.3, at = 2, cex = 1, adj = 0)
mtext(text = "Purifying", side = 2, line = 2.3, at = -2, cex = 1, adj = 1)
x.lim
arrows(x0 = -2, y0 = 2, x1 = -2, y1 = 17, xpd = T, length = .075)
arrows(x0 = -2, y0 = -2, x1 = -2, y1 = -19, xpd = T, length = .075)
dev.off()
# ----