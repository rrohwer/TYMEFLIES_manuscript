# RRR
# DOC boxplot
# fig 4i

library(data.table)

doc <- fread(input = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/dissolved_organic_carbon.csv")

# ---- make plot ----

doc <- doc[!is.na(Mean.DOC.mg.L)]
doc[ , color.2012 := adjustcolor("grey40", .7)]
doc[year4 == 2012 , color.2012 := "red2"]
doc[ , cex.2012 := .8]
doc[year4 == 2012 , cex.2012 := 1.4]


y.lim <- c(min(doc$Mean.DOC.mg.L, na.rm = T), max(doc$Mean.DOC.mg.L, na.rm = T))
y.ticks <- seq.int(1,8,2)



plot(x = c(.6,1.4), y = y.lim, type = "n", ann = F, axes = F)
points(x = runif(n = nrow(doc[year4 != 2012]), min = .7, max = 1.3), y = doc[year4 != 2012, Mean.DOC.mg.L], col = doc[year4 != 2012, color.2012], cex = doc[year4 != 2012, cex.2012], pch = 16)
boxplot(x = doc[ , Mean.DOC.mg.L], pch = NA, lty = 1, col = adjustcolor("white",0), ylim = y.lim, axes = F, add = T, at = 1)
points(x = runif(n = nrow(doc[year4 == 2012]), min = .9, max = 1.1), y = doc[year4 == 2012, Mean.DOC.mg.L], col = doc[year4 == 2012, color.2012], cex = doc[year4 == 2012, cex.2012], pch = 16)
box()
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.075)
axis(side = 2, at = y.ticks, labels = y.ticks, lwd = 0, line = -.5, las = 2)
mtext(text = "DOC (mg/L)", side = 2, line = 1.5)
