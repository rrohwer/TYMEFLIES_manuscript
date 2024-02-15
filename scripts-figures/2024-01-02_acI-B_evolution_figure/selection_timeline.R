# RRR
# data for the acI-B figure: ME2011-09-21_3300043464_group3_bin69
# Fig 5d - selection over time

library(data.table)
library(lubridate)

sel <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/amt_selection-ME2011-09-21_3300043464_group3_bin69.csv", colClasses = c("date" = "character"))

index.change <- which(sel$date == "2012-08-03")

sel[ ,date := parse_date_time(date, "ymd")]

change.date <- sel$date[index.change]

fill.under.lines <- function(X, Y, YAxisMin, Color, xpd = F){
  poly.x <- c(min(X), X, max(X))
  poly.y <- c(YAxisMin, Y, YAxisMin )
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

x.ticks <- parse_date_time(paste(seq(2000,2020,2), 1, 1), "ymd")
x.labs <- year(x.ticks)
x.lab.locs <- x.ticks
x.lab.locs[x.labs == 2020] <- x.lab.locs[x.labs == 2020] %m+% -months(6)

summary(sel$Npos)
y.ticks <- seq(0,17,4)

plot(x = sel[ ,date], y = sel[ ,Npos], type = "n", ann = F, axes = F, xlim = c(min(x.ticks), max(x.ticks)), xaxs = "i")
fill.under.lines(X = sel[ ,date], Y = sel[ ,Npos], YAxisMin = 0, Color = adjustcolor("black",.2))
abline(v = change.date, col = "black", xpd = F, lwd = 1)
lines(x = sel[ ,date], y = sel[ ,Npos], lwd = .7, col = adjustcolor("black",.5))
points(x = sel[ ,date], y = sel[ ,Npos], pch = 16, col = adjustcolor(sel$color.year.step,.7))
box()
axis(side = 1, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03, xpd = NA)
axis(side = 1, at = x.lab.locs, labels = x.labs, lwd = 0, line = -.55, xpd = NA)
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 2, at = y.ticks, labels = T, lwd = 0, las = 2, line = -.3)
mtext(text = "Positively Selected\nGenes (count)", side = 2, line = 2.1, at = 8)

