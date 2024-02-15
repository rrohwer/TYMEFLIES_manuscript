# RRR
# show off good completeness and very very low contamination

library(data.table)

checkm <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv")
completeness <- checkm$completeness
contamination <- checkm$contamination

plot(x = 1, y = 1, xlim = c(.75,1.25), ylim = c(50,100), type = "n", ann = F, axes = F)
points(x = runif(n = length(completeness), min = .75, max = 1.25), y = completeness, pch = 16, col = adjustcolor("slategray",.2), cex = .4)
boxplot(x = completeness, lty = 1, col = adjustcolor("white",0), add = T, at = 1, lwd = 2, axes = F, ann = F)
box()
axis(side = 2, at = seq(50,100,25), labels = F, lwd = 0, lwd.ticks = 1, tck = -.075, line = 0)
axis(side = 2, at = seq(50,100,25), labels = T, lwd = 0, las = 2, line = -.5)
mtext(text = "Completeness (%)", side = 2, line = 1.75)
mtext(text = "Total Genomes: 2,855", side = 1, line = .4, at = .87, adj = 0)
