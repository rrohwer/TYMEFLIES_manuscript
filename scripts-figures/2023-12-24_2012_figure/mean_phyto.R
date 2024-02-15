# RRR
# phyto boxplot
# fig 4h

library(data.table)

phyto <- fread(input = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/phyto_mean_biomass.csv")

# ---- make plot ----

phyto <- phyto[!is.na(Mean.Biomass.mg.L)]
phyto[ , color.2012 := adjustcolor("grey40", .7)]
phyto[year4 == 2012 , color.2012 := "red2"]
phyto[ , cex.2012 := .8]
phyto[year4 == 2012 , cex.2012 := 1.4]


y.lim <- c(min(phyto$Mean.Biomass.mg.L, na.rm = T), max(phyto$Mean.Biomass.mg.L, na.rm = T))
y.ticks <- seq.int(0,10,3)


plot(x = c(.6,1.4), y = y.lim, type = "n", ann = F, axes = F)
points(x = runif(n = nrow(phyto[year4 != 2012]), min = .7, max = 1.3), y = phyto[year4 != 2012, Mean.Biomass.mg.L], col = phyto[year4 != 2012, color.2012], cex = phyto[year4 != 2012, cex.2012], pch = 16)
boxplot(x = phyto[ , Mean.Biomass.mg.L], pch = NA, lty = 1, col = adjustcolor("white",0), ylim = y.lim, axes = F, add = T, at = 1)
points(x = runif(n = nrow(phyto[year4 == 2012]), min = .9, max = 1.1), y = phyto[year4 == 2012, Mean.Biomass.mg.L], col = phyto[year4 == 2012, color.2012], cex = phyto[year4 == 2012, cex.2012], pch = 16)
box()
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.075)
axis(side = 2, at = y.ticks, labels = y.ticks, lwd = 0, line = -.5, las = 2)
mtext(text = "Phytoplankton (mg/L)", side = 2, line = 1.5)

