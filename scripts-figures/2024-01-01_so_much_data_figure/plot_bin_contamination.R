# RRR
# show off good completeness and very very low contamination

library(data.table)

checkm <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv")
contamination <- checkm$contamination
contamination <- checkm$contamination

plot(x = 1, y = 1, xlim = c(.75,1.25), ylim = c(0,10), type = "n", ann = F, axes = F)
points(x = runif(n = length(contamination), min = .75, max = 1.25), y = contamination, pch = 16, col = adjustcolor("slategray",.2), cex = .4)
boxplot(x = contamination, lty = 1, col = adjustcolor("white",0), add = T, at = 1, lwd = 2, axes = F, ann = F, pch = NA)
box()
axis(side = 2, at = seq(0,10,5), labels = F, lwd = 0, lwd.ticks = 1, tck = -.075, line = 0)
axis(side = 2, at = seq(0,10,5), labels = T, lwd = 0, las = 2, line = -.5)
mtext(text = "Contamination (%)", side = 2, line = 1.3)

