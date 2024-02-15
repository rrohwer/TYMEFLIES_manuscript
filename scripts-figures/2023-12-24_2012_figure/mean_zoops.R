# RRR
# mean non-pred zoop biomass boxplot
# for fig 4 panel d

library(data.table)

zoop <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/zoop_mean_annual_biomass.csv")

# ---- make plot ----

zoop <- zoop[!is.na(Mean.Biomass.mg.L)]
zoop[ , color.2012 := adjustcolor("grey40", .7)]
zoop[year4 == 2012 , color.2012 := "red2"]
zoop[ , cex.2012 := .8]
zoop[year4 == 2012 , cex.2012 := 1.4]


y.lim <- c(min(zoop$Mean.Biomass.mg.L, na.rm = T), max(zoop$Mean.Biomass.mg.L + .5, na.rm = T))
y.ticks <- seq.int(0,16,4)


plot(x = c(.6,1.4), y = y.lim, type = "n", ann = F, axes = F)
points(x = runif(n = nrow(zoop[year4 != 2012]), min = .7, max = 1.3), y = zoop[year4 != 2012, Mean.Biomass.mg.L], col = zoop[year4 != 2012, color.2012], cex = zoop[year4 != 2012, cex.2012], pch = 16)
boxplot(x = zoop[ , Mean.Biomass.mg.L], pch = NA, lty = 1, col = adjustcolor("white",0), ylim = y.lim, axes = F, add = TRUE, at = 1)
points(x = runif(n = nrow(zoop[year4 == 2012]), min = .9, max = 1.1), y = zoop[year4 == 2012, Mean.Biomass.mg.L], col = zoop[year4 == 2012, color.2012], cex = zoop[year4 == 2012, cex.2012], pch = 16)
box()
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.075)
axis(side = 2, at = y.ticks, labels = y.ticks, lwd = 0, line = -.5, las = 2)
mtext(text = "Zooplankton (mg/L)", side = 2, line = 2)
