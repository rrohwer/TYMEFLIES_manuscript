# RRR
# fig 5e

library(data.table)

nmds <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/nmds_coords-ME2011-09-21_3300043464_group3_bin69.csv", colClasses = c("date" = "character"))

background.index <- which(nmds$Year != 2011 & nmds$Year != 2012 & nmds$Year != 2013 & nmds$Year != 2014)
forground.index  <- which(nmds$Year == 2011 | nmds$Year == 2013 | nmds$Year == 2014)
index.2012 <- which(nmds$Year == 2012)

x.lim <- c(-40,60) # range of 100
y.lim <- c(-25,25) # range of 50
x.ticks <- seq(-40,55,12.5)
y.ticks <- seq(-25,25,12.5)

plot(x = nmds$x, y = nmds$y, xlim = x.lim, ylim = y.lim, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .5)
mtext(text = "SNV Profiles, NMDS Axis 2", side = 2, line = .5)
axis(side = 1, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = .03)
axis(side = 3, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = .03)
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = .03)
axis(side = 4, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = .03)
points(x = nmds[background.index, x], y = nmds[background.index, y], pch = 16, cex = 1.3, col = adjustcolor(nmds[background.index, color.year.step], alpha.f = .7))
points(x = nmds[forground.index, x], y = nmds[forground.index, y], pch = 16, cex = 1.3, col = adjustcolor(nmds[forground.index, color.year.step], alpha.f = .7))
points(x = nmds[index.2012, x], y = nmds[index.2012, y], pch = 16, cex = 1.3, col = adjustcolor(nmds[index.2012, color.year.step], alpha.f = .9))

# text(x = nmds[index.2012, x], y = nmds[index.2012, y], labels = nmds[index.2012,date])

change.index <- which(nmds$date == "2012-08-03")
points(x = nmds[change.index, x], y = nmds[change.index, y], pch = 21, cex = 3)
# segments(x0 = nmds[change.index, x], x1 = nmds[change.index, x] + 3, y0 = nmds[change.index, y], y1 = nmds[change.index, y] - 4)
text(x = nmds[change.index, x] + 1.5, y = nmds[change.index, y] - 4, labels = "2012-08-03", adj = 0, srt = 0)


