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

# acI.C <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2016-07-20_3300033996_group7_bin32_gene_info_combined_MK.tsv.gz")
# genes <- acI.C
# my.genome <- "ME2016-07-20_3300033996_group7_bin32"
# 
# acI.A <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-04_3300044729_group3_bin142_gene_info_combined_MK.tsv.gz")

# ---- parse data ----

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,sample := sub(pattern = "\\.IS_gene_info\\.tsv", replacement = "", x = sample)]
genes[ ,sample := sub(pattern = ".gz", replacement = "", x = sample)]

pos <- genes[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman), .(.N), by = .(sample, date, season, year)]
neg <- genes[MK.p.val <= .05 & mcdonald.kreitman > 1 & is.finite(mcdonald.kreitman), .(.N), by = .(sample, date, season, year)] # THIS IS A BUG MAKING THESE PLOTS INCORRECT- pvalue should be less than .05!!!
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


# # for acI-C ----
# pdf(file = "figures/2023-08-05_ESA_plots/gene_selection_acI-C_seasons.pdf", width = 6.5, height = 3)
# bar.spots <- barplot(height = season.sum$Npos.per.sample, beside = T, plot = F, space = .75)
# x.lim <- c(min(bar.spots), max(bar.spots))
# y.lim <- c(-max(season.sum$Nneg.per.sample, na.rm = T), max(season.sum$Npos.per.sample, na.rm = T))
# par(mar = c(2.5,6,2,1))
# plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
# barplot(height = season.sum$Npos.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = "#b89996", add = T, border = NA, space = .75)
# barplot(height = -season.sum$Nneg.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = "#7592b2", add = T, border = NA, space = .75)
# segments(x0 = .3, x1 = 11, y0 = 0, y1 = 0, xpd = T)
# text(x = bar.spots, y = -45, labels = season.sum$label, cex = 1, xpd = NA)
# y.lim
# ticks <- axis(side = 2, line = 1, labels = F, tck = -.02, at = seq.int(-40,10,10), xpd = T)
# axis(side = 2, line = .75, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = 1, xpd = T)
# mtext(text = "Genes Under Selection (ave / sample)", side = 2, line = 5, cex = 1)
# mtext(text = "Positive", side = 2, line = 3.5, at = 2, cex = 1, adj = 0)
# mtext(text = "Purifying", side = 2, line = 3.5, at = -2, cex = 1, adj = 1)
# x.lim
# arrows(x0 = -.3, y0 = 2, x1 = -.3, y1 = 17, xpd = T, length = .075)
# arrows(x0 = -.3, y0 = -2, x1 = -.3, y1 = -19, xpd = T, length = .075)
# dev.off()
# # ----
# for acI-B ----
pdf(file = "figures/2023-08-10_UCSB_plots/gene_selection_acI-B_long-term-positive_only_for_UGA.pdf", width = 3.1, height = 1.36)
bar.spots <- barplot(height = year.sum$Npos.per.sample, beside = T, plot = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(year.sum$Nneg.per.sample, na.rm = T), max(year.sum$Npos.per.sample, na.rm = T))
x.ticks <- parse_date_time(x = c(2000,2012,2019), orders = "y")
par(mar = c(.9,3.4,.38,.1))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F, yaxs = "i")
barplot(height = year.sum$Npos.per.sample, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = adjustcolor("#765a5a", alpha.f = 1, red.f = 1.5), add = T, border = NA)
axis(side = 1, at = c(min(bar.spots)-1, max(bar.spots)+.5), lwd = 1, lwd.ticks = 0, labels = F, line = -.305, xpd = NA)
axis(side = 1, at = bar.spots[c(1,13,20)], labels = F, tck = -.05, lwd = 0, lwd.ticks = 1, line = -.305)
axis(side = 1, at = bar.spots[c(1)], labels = year(x.ticks)[1], lwd = 0, las = 1, line = -1.05, hadj = .25)
axis(side = 1, at = bar.spots[c(13)], labels = year(x.ticks)[2], lwd = 0, las = 1, line = -1.05, hadj = .5)
axis(side = 1, at = bar.spots[c(20)], labels = year(x.ticks)[3], lwd = 0, las = 1, line = -1.05, hadj = .75)

ticks <- axis(side = 2, line = 0, labels = F, tck = -.05, at = c(0,4,8), xpd = T)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = 1, xpd = T)
mtext(text = "Genes Under\nSelection", side = 2, line = 1.4, cex = 1)

bar.spots
rect(xleft = bar.spots[2], xright = bar.spots[3], ybottom = -26, ytop = -18, col = adjustcolor("#765a5a", alpha.f = 1, red.f = 1.5), border = NA)
text(x = bar.spots[4] - .5, y = -22, labels = "Positive", adj = 0)
rect(xleft = bar.spots[2], xright = bar.spots[3], ybottom = -38, ytop = -30, col = adjustcolor("#a9bfdb", alpha.f = 1), border = NA)
text(x = bar.spots[4] - .5, y = -34, labels = "Purifying", adj = 0)

dev.off()
# ----