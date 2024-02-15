# RRR

abund <- readRDS("figures/2023-12-24_paper_figures/Rohwer_Figure_1_data/16S_rank_abund.rds")



bar.spots <- barplot(height = abund, horiz = T, axes = F, ann = F, names.arg = rep("",ncol(abund)), space = .5, col = c("slategray","orange2"))
text(x = -1, y = bar.spots, labels = colnames(abund), adj = 1, xpd = NA)
axis(side = 3, at = seq(0,28,14), labels = F, tck = -.03, line = 0, xpd = NA)
axis(side = 3, at = seq(0,28,14), labels = T, lwd = 0, line = -.65, xpd = NA)
mtext(text = "16S Abundance (%)", side = 3, line = 1.5, xpd = NA)

# For text: mean abundance of Frankiales is 22 %
abund[2,7]
