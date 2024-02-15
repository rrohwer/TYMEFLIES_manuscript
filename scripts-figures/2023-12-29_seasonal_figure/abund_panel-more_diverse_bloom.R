# RRR 
# abundance panel of the less diverse blooms example plot
# fig 2c, top
# example genomes: 
# less diverse bloom/anicorrelated: ME2011-09-04_3300044729_group3_bin142
# more diverse bloom/correlated: ME2012-08-31_3300044613_group4_bin150

library(data.table)
library(lubridate)

abund <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/more_diverse_bloom-ME2012-08-31_3300044613_group4_bin150.csv")
ave <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/more_diverse_bloom_1mo_moving_average-ME2012-08-31_3300044613_group4_bin150.csv")

fill.btwn.lines <- function(X, Y1, Y2, Color, xpd = F){
  index <- is.na(X) | is.na(Y1) | is.na(Y2)
  X <- X[!index]
  Y1 <- Y1[!index]
  Y2 <- Y2[!index]
  poly.x <- c(X, X[length(X):1])
  poly.y <- c(Y1, Y2[length(Y2):1])
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

x.ticks <- parse_date_time(x = seq(2,12,2), orders = "m")
x.labs <- month(x.ticks, abbr = T, label = T)
x.ticks <- yday(x.ticks)
x.lab.locs <- x.ticks 

max(abund$abund.perc)
y.lim <- c(0,1.25)
y.ticks <- seq(0,1.25,.5)

plot(x = abund$yday, y = abund$abund.perc, type = "n", ann = F, axes = F, xlim = c(0,365), ylim = y.lim, xaxs = "i")
fill.btwn.lines(X = ave[num.obs.abund > 0, center], Y1 = ave[num.obs.abund > 0, ave.abund + sd.abund], Y2 = ave[num.obs.abund > 0, ave.abund - sd.abund], Color = adjustcolor("grey40", .4))
for (yr in unique(abund$year)){
  points(x = abund[year == yr, yday], y = abund[year == yr, abund.perc], pch = 16, cex = .4, col = adjustcolor("royalblue",1))
  lines(x = abund[year == yr, yday], y = abund[year == yr, abund.perc], lwd = .5, col = adjustcolor("royalblue",.5))
}
lines(x = ave[num.obs.abund > 0, center], y = ave[num.obs.abund > 0, ave.abund], lwd = 3, col = "black")
box()
axis(side = 1, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
# axis(side = 1, at = x.lab.locs, labels = x.labs, lwd = 0, line = -.5)
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 2, at = y.ticks, labels = y.ticks, lwd = 0, line = -.5, las = 2)
mtext(text = "Abundance (%)", side = 2, line = 2.75, xpd = NA)







