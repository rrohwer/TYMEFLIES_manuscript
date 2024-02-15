# RRR
# SRP boxplot
# fig 4f

library(data.table)

srp <- fread(input = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/soluble_reactive_phosphorus.csv")

# ---- make plot ----

srp <- srp[!is.na(Mean.SRP.mg.L)]
srp[ , color.2012 := adjustcolor("grey40", .7)]
srp[year4 == 2012 , color.2012 := "red2"]
srp[ , cex.2012 := .8]
srp[year4 == 2012 , cex.2012 := 1.4]


y.lim <- c(min(srp$Mean.SRP.mg.L, na.rm = T), max(srp$Mean.SRP.mg.L, na.rm = T))
y.ticks <- seq.int(0.02,0.09,.02)


plot(x = c(.6,1.4), y = y.lim, type = "n", ann = F, axes = F)
points(x = runif(n = nrow(srp[year4 != 2012]), min = .7, max = 1.3), y = srp[year4 != 2012, Mean.SRP.mg.L], col = srp[year4 != 2012, color.2012], cex = srp[year4 != 2012, cex.2012], pch = 16)
boxplot(x = srp[ , Mean.SRP.mg.L], pch = NA, lty = 1, col = adjustcolor("white",0), ylim = y.lim, axes = F, add = T, at = 1)
points(x = runif(n = nrow(srp[year4 == 2012]), min = .9, max = 1.1), y = srp[year4 == 2012, Mean.SRP.mg.L], col = srp[year4 == 2012, color.2012], cex = srp[year4 == 2012, cex.2012], pch = 16)
box()
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.075)
axis(side = 2, at = y.ticks, labels = y.ticks, lwd = 0, line = -.5, las = 2)
mtext(text = "SRP (mg/L)", side = 2, line = 2.5)
