# RRR

abund <- readRDS("figures/2023-12-24_paper_figures/Rohwer_Figure_1_data/MAG_rank_abund.rds")



bar.spots <- barplot(height = abund, horiz = T, axes = F, ann = F, names.arg = rep("",ncol(abund)), space = .5, col = c("slategray","orange2"))

axis(side = 3, at = seq(0,16,8), labels = F, tck = -.03, line = 0, xpd = NA)
axis(side = 3, at = seq(0,16,8), labels = T, lwd = 0, line = -.65, xpd = NA)
mtext(text = "MAG Abundance (%)", side = 3, line = 1.5, xpd = NA)

text(x = 7, y = bar.spots[3] - (bar.spots[3] - bar.spots[2]) / 2, labels = "Nanopelagicales\nOrder", adj = 0, xpd = NA)
rect(xleft = 5, xright = 6.25, ybottom = bar.spots[2], ytop = bar.spots[3], col = "orange2")

# For text: mean abundance of Frankiales is 10 %
abund[2,7]
