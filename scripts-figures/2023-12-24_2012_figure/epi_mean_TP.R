# RRR
# TP boxplot
# fig 4f

library(data.table)

tp <- fread(input = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/total_phosphorus.csv")

# ---- make plot ----

tp <- tp[!is.na(Mean.TP.mg.L)]
tp[ , color.2012 := adjustcolor("grey40", .7)]
tp[year4 == 2012 , color.2012 := "red2"]
tp[ , cex.2012 := .8]
tp[year4 == 2012 , cex.2012 := 1.4]


y.lim <- c(min(tp$Mean.TP.mg.L, na.rm = T) - .002, max(tp$Mean.TP.mg.L, na.rm = T))
y.ticks <- seq.int(0.05,0.12,.02)


plot(x = c(.6,1.4), y = y.lim, type = "n", ann = F, axes = F)
points(x = runif(n = nrow(tp[year4 != 2012]), min = .7, max = 1.3), y = tp[year4 != 2012, Mean.TP.mg.L], col = tp[year4 != 2012, color.2012], cex = tp[year4 != 2012, cex.2012], pch = 16)
boxplot(x = tp[ , Mean.TP.mg.L], pch = NA, lty = 1, col = adjustcolor("white",0), ylim = y.lim, axes = F, add = T, at = 1)
points(x = runif(n = nrow(tp[year4 == 2012]), min = .9, max = 1.1), y = tp[year4 == 2012, Mean.TP.mg.L], col = tp[year4 == 2012, color.2012], cex = tp[year4 == 2012, cex.2012], pch = 16)
box()
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.075)
axis(side = 2, at = y.ticks, labels = y.ticks, lwd = 0, line = -.5, las = 2)
mtext(text = "TP (mg/L)", side = 2, line = 2.5)

