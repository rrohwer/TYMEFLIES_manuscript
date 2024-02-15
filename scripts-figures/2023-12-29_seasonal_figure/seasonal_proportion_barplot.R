# RRR
# what proportion of all mags had seasonality in pi and abundance?
# barplot split by taxonomy

seas <- readRDS("figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/seasonality_summary_data.rds")


y.lim <- max(colSums(seas))
bar.spots <- barplot(height = seas, legend.text = FALSE, col = c("grey80","#e79493","#9865c0","#a9adff"), axes = F, ann = F, names.arg = rep("",ncol(seas)))
text(x = bar.spots, y = -(y.lim / 20), labels = colnames(seas), xpd = NA, srt = 30, adj = 1)
axis(side = 2, at = c(0,y.lim), labels = F, lwd = 1, lwd.ticks = 0, line = -.3)
axis(side = 2, at = seq(0,100,25), labels = F, line = -.3, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 2, at = seq(0,100,25), labels = T, line = -.55, lwd = 0, lwd.ticks = 0, las = 2)
mtext(text = "Percent Genomes", side = 2, line = 2, outer = F, at = 60)
# mtext(text = "Total Genomes: 1,474", side = 2, line = 2, outer = F, at = 60) # not fitting well, put in caption

# # side legend
# y.legend.locs <- 35 - (y.lim / 13) * (1:5)
# bar.spacing <- bar.spots[2] - bar.spots[1]
# rect(xleft = 0 - bar.spacing * 3.7, xright = 0 - bar.spacing * 3.7 + bar.spacing / 3, ytop = y.legend.locs[-1] - (y.lim / 40), ybottom = y.legend.locs[-1] + (y.lim / 40), xpd = NA, col = c("#a9adff","#9865c0","#e79493","grey80"))
# text(x = 0 - bar.spacing * 3.7 + bar.spacing / 3 + bar.spacing / 10, y = y.legend.locs[-1], labels = c("Diversity","Both","Abund.","Neither"), xpd = NA, adj = 0)
# text(x = 0 - bar.spacing * 3.7, y = y.legend.locs[1] + 1, labels = "Seasonal:", xpd = NA, adj = 0)

# Top legend (flat)
y.legend.locs <- 112
half.bar.height <- 3
bar.spacing <- bar.spots[2] - bar.spots[1]
bar.width <- bar.spacing / 3.25
x.legend.locs <- (0 - bar.spacing * 2.25) + seq(0, bar.spots[8] + bar.width * 3, along.with = 1:5)
text.gap <- bar.spacing / 7
x.legend.locs <- x.legend.locs + c(text.gap * 6, 0, text.gap*1, text.gap * 5, + text.gap * .75 )

rect(xleft = x.legend.locs[-1], xright = x.legend.locs[-1] + bar.width, ytop = y.legend.locs - 2 + half.bar.height, ybottom = y.legend.locs - 2 - half.bar.height, xpd = NA, col = c("#a9adff","#e79493","#9865c0","grey80"))
text(x = x.legend.locs[-1] + bar.width + text.gap, y = y.legend.locs - 2, labels = c("Diversity","Abundance","Both","Neither"), xpd = NA, adj = 0)
text(x = x.legend.locs[1] + bar.width, y = y.legend.locs + 8, labels = "Seasonality in:", xpd = NA, adj = 0)
