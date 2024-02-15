# RRR
# what proportion of all mags with seasonality had more/less diverse bloom patterns?
# barplot split by taxonomy

seas <- readRDS("figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/bloom_diversity_summary_data.rds")

# # check
# barplot(seas,legend.text = T, col = c("grey80",adjustcolor("green4",.7),adjustcolor("gold2",.7)))

y.lim <- max(colSums(seas))
bar.spots <- barplot(height = seas, legend.text = FALSE, col = c("grey80",adjustcolor("green4",.7),adjustcolor("gold2",.7)), axes = F, ann = F, names.arg = rep("",ncol(seas)))
text(x = bar.spots, y = -(y.lim / 20), labels = colnames(seas), xpd = NA, srt = 30, adj = 1)
axis(side = 2, at = c(0,y.lim), labels = F, lwd = 1, lwd.ticks = 0, line = -.3)
axis(side = 2, at = seq(0,100,25), labels = F, line = -.3, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 2, at = seq(0,100,25), labels = T, line = -.55, lwd = 0, lwd.ticks = 0, las = 2)
mtext(text = "Percent Genomes", side = 2, line = 2, outer = F, at = 60)

# # side legend
# y.legend.locs <- 35 - (y.lim / 13) * (1:5)
# bar.spacing <- bar.spots[2] - bar.spots[1]
# rect(xleft = 0 - bar.spacing * 3.7, xright = 0 - bar.spacing * 3.7 + bar.spacing / 3, ytop = y.legend.locs[-1] - (y.lim / 40), ybottom = y.legend.locs[-1] + (y.lim / 40), xpd = NA, col = c("#a9adff","#9865c0","#e79493","grey80"))
# text(x = 0 - bar.spacing * 3.7 + bar.spacing / 3 + bar.spacing / 10, y = y.legend.locs[-1], labels = c("Diversity","Both","Abund.","Neither"), xpd = NA, adj = 0)
# text(x = 0 - bar.spacing * 3.7, y = y.legend.locs[1] + 1, labels = "Seasonal:", xpd = NA, adj = 0)

# Top legend (flat)
y.legend.locs <- 114
half.bar.height <- 3
bar.spacing <- bar.spots[2] - bar.spots[1]
bar.width <- bar.spacing / 3
x.legend.locs <- (1 - bar.spacing * 1.75) + seq(0, bar.spots[6] + bar.width * 3, along.with = 1:4)
text.gap <- bar.spacing / 5
x.legend.locs <- x.legend.locs + c(-text.gap * 3, 0,0,0)

rect(xleft = x.legend.locs[-1], xright = x.legend.locs[-1] + bar.width, ytop = y.legend.locs + half.bar.height, ybottom = y.legend.locs - half.bar.height, xpd = NA, col = c(adjustcolor("gold2",.7), adjustcolor("green4",.7), "grey80"))
text(x = x.legend.locs[-1] + bar.width + text.gap, y = y.legend.locs, labels = c("Less\nDiverse", "More\nDiverse", "Not\nCorrelated"), xpd = NA, adj = 0)
text(x = x.legend.locs[1] + bar.width, y = y.legend.locs, labels = "Blooms are:", xpd = NA, adj = 0)
