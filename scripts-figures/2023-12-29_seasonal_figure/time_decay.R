# RRR
# make the time decay figure
# fig 2A

library(data.table)

dist.each <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/time_decay_distances-ME2017-06-13_3300043469_group7_bin14.csv.gz")
moving.ave <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/time_decay_6mo_moving_average-ME2017-06-13_3300043469_group7_bin14.csv.gz")


fill.btwn.lines <- function(X, Y1, Y2, Color, xpd = F){
  index <- is.na(X) | is.na(Y1) | is.na(Y2)
  X <- X[!index]
  Y1 <- Y1[!index]
  Y2 <- Y2[!index]
  poly.x <- c(X, X[length(X):1])
  poly.y <- c(Y1, Y2[length(Y2):1])
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}


max(dist.each$dist) # 59.54591
min(dist.each$dist) # 10.51186

plot(dist ~ time.approx.year, data = dist.each, type = "n", ann = F, axes = F, xaxs = "i", xlim = c(-.1,19.25), yaxs = "i", ylim = c(10,60))
# fill.btwn.lines(X = moving.ave$time.approx.year, Y1 = moving.ave$ave + moving.ave$sd, Y2 = moving.ave$ave - moving.ave$sd, Color = adjustcolor("black",.2))
points(dist ~ time.approx.year, data = dist.each, pch = 16, col = adjustcolor("royalblue",.3), cex = .3)
lines(dist ~ time.approx.year, data = moving.ave[num.obs >= 100], lwd = 3, col = "black")
box()
axis(side = 1, at = seq(0,19,2), labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 1, at = seq(0,19,2), lwd = 0, line = -.5)
axis(side = 2, at = seq(10,60,20), labels = FALSE, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 2, at = seq(10,60,20), las = 2, lwd = 0, line = -.35)
# mtext(text = "Euclidean Distance", side = 2, line = 3)
mtext(text = "Distance Between SNV Profiles", side = 2, line = 2)
mtext(text = "Years between timepoints", side = 1, line = 1.65)

points(x = 14.25, y = 56.5, pch = 19, col = adjustcolor("royalblue",.7), cex = .75)
text(x = 14.75, y = 56.5, labels = "Pairwise comparison", adj = 0)
segments(x0 = 14.10, x1 = 14.40, y0 = 52, y1 = 52, lwd = 3, col = "black")
text(x = 14.75, y = 52, labels = "Moving average", adj = 0)
# rect(xleft = 13.85, xright = 19.25, ybottom = 48, ytop = 60)






