# RRR
# ice duration boxplot
# fig 4c

library(data.table)

ice <- fread(input = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/ice_durations.csv")

# ---- make plot ----

ice <- ice[!is.na(ice_duration)]
ice[ , color.2012 := adjustcolor("grey40", .7)]
ice[summer.year == 2012 , color.2012 := "red2"]
ice[ , cex.2012 := .8]
ice[summer.year == 2012 , cex.2012 := 1.4]


y.lim <- c(min(ice$ice_duration, na.rm = T), max(ice$ice_duration, na.rm = T))
y.ticks <- seq.int(20,160,35)


plot(x = c(.6,1.4), y = y.lim, type = "n", ann = F, axes = F)
points(x = runif(n = nrow(ice[summer.year != 2012]), min = .7, max = 1.3), y = ice[summer.year != 2012, ice_duration], col = ice[summer.year != 2012, color.2012], cex = ice[summer.year != 2012, cex.2012], pch = 16)
boxplot(x = ice[ , ice_duration], pch = NA, lty = 1, col = adjustcolor("white",0), ylim = y.lim, axes = F, add = TRUE, at = 1)
points(x = runif(n = nrow(ice[summer.year == 2012]), min = .9, max = 1.1), y = ice[summer.year == 2012, ice_duration], col = ice[summer.year == 2012, color.2012], cex = ice[summer.year == 2012, cex.2012], pch = 16)
box()
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.075)
axis(side = 2, at = y.ticks, labels = y.ticks, lwd = 0, line = -.5, las = 2)
mtext(text = "Ice Duration (days)", side = 2, line = 2.5)
