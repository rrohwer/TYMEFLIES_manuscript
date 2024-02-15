# RRR
# make fig 4a, the abrupt SNV changes plot
# try it as a beeswarm plot instead of a jittered plot

# beeswarm doesn't have a way to adjust the point cex as part of its default plot
# also, it changes the order of the labels to be alphabetical, but I need them arranged by upper taxonomy
# so have to get coords for each genus separately and then piece together

library(data.table)
library(pBrackets)

snvs <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/abrupt_genomic_changes.csv")

# ---- save plotting info ----

x.locs <- unique(snvs$label.loc)
x.labs <- unique(snvs$group.ID)
x.col <- unique(snvs$phylum.color)

plot(x = c(min(x.locs)-.8,max(x.locs)+.5), y = c(2000,2020), type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
abline(v = x.locs, col = adjustcolor("black",.2))
# abline(h = 2012:2013, col = adjustcolor("black",.2))
box()
# axis(side = 1, at = x.locs, labels = F, lwd = 0, lwd.ticks = 1)
text(x = x.locs, y = 1999.25, labels = x.labs, adj = 1, xpd = NA, srt = 90, col = x.col)
axis(side = 2, at = seq.int(2000,2019,2), labels = F, lwd = 0, lwd.ticks = 1, tck = -.02)
axis(side = 2, at = seq.int(2000,2019,2), labels = T, lwd = 0, line = -.5, las = 2)
points(x = snvs$label.loc + snvs$x.shift, y = snvs$change.date, pch = 16, col = "white", cex = snvs$cex) 
points(x = snvs$label.loc + snvs$x.shift, y = snvs$change.date, pch = 16, col = adjustcolor(snvs$color,.7), cex = snvs$cex) 
text(x = 2, y = 1991.25, labels = "Actinobacteriota", xpd = NA)
brackets(x2 = .7, x1 = 3.3, y1 = 1992.75, y2 = 1992.75, xpd = NA, h = .5, type = 4, ticks = NA)

rect(xleft = 5, xright = max(x.locs)+.5, ybottom = 2016.75, ytop = 2020, col = "white")
points(x = c(5.3,5.3), y = c(2019,2017.75), pch = 16, col = adjustcolor(unique(snvs$color),.7), cex = 1.5)
text(x = c(5.6,5.6), y = c(2019,2017.75), labels = c("Step change","Disturbance/resilience"), adj = 0)
