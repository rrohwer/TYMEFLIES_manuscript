# RRR 
# diversity panel of the less diverse blooms example plot
# fig 2c, bottom
# example genomes: 
# less diverse bloom/anicorrelated: ME2011-09-04_3300044729_group3_bin142
# more diverse bloom/correlated: ME2012-08-31_3300044613_group4_bin150

library(data.table)
library(lubridate)

div <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/more_diverse_bloom-ME2012-08-31_3300044613_group4_bin150.csv")
ave <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/more_diverse_bloom_1mo_moving_average-ME2012-08-31_3300044613_group4_bin150.csv")

div <- div[!is.na(nucl_diversity)]
ave <- ave[!is.na(ave.div)]

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

max(div$nucl_diversity)
min(div$nucl_diversity)
y.lim <- c(.009,.022)
y.ticks <- seq(0,.03,.005)

plot(x = div$yday, y = div$nucl_diversity, type = "n", ann = F, axes = F, xlim = c(0,365), ylim = y.lim, xaxs = "i")
fill.btwn.lines(X = ave[num.obs.div > 0, center], Y1 = ave[num.obs.div > 0, ave.div + sd.div], Y2 = ave[num.obs.div > 0, ave.div - sd.div], Color = adjustcolor("grey40", .4))
for (yr in unique(div$year)){
  points(x = div[year == yr, yday], y = div[year == yr, nucl_diversity], pch = 16, cex = .4, col = adjustcolor("royalblue",1))
  lines(x = div[year == yr, yday], y = div[year == yr, nucl_diversity], lwd = .5, col = adjustcolor("royalblue",.5))
}
lines(x = ave[num.obs.div > 0, center], y = ave[num.obs.div > 0, ave.div], lwd = 3, col = "black")
box()
axis(side = 1, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 1, at = x.lab.locs, labels = x.labs, lwd = 0, line = -.65, xpd = T)
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 2, at = y.ticks, labels = y.ticks, lwd = 0, line = -.5, las = 2)
mtext(text = expression(paste("Diversity (",pi,")")), side = 2, line = 2.75, xpd = NA)








