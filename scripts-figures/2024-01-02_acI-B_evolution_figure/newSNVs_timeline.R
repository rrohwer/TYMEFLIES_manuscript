# RRR
# data for the acI-B figure: ME2011-09-21_3300043464_group3_bin69
# Fig 5c - new SNVs over time

library(data.table)
library(lubridate)

snv <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/new_snvs_ME2011-09-21_3300043464_group3_bin69.csv", colClasses = c("Date" = "character"))

index.change <- which(snv$Date == "2012-08-03")

snv[ ,Date := parse_date_time(Date, "ymd")]

change.date <- snv$Date[index.change]

fill.under.lines <- function(X, Y, YAxisMin, Color, xpd = F){
  poly.x <- c(min(X), X, max(X))
  poly.y <- c(YAxisMin, Y, YAxisMin )
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

x.ticks <- parse_date_time(paste(seq(2000,2020,2), 1, 1), "ymd")
# x.labs <- year(x.ticks)
# x.lab.locs <- x.ticks
# x.lab.locs[x.labs == 2020] <- x.lab.locs[x.labs == 2020] %m+% -months(0)

# truncate axis so initial high points are excluded
summary(snv$New)
y.ticks <- seq(0,1500,500)

plot(x = snv[ ,Date], y = snv[ ,New], type = "n", ann = F, axes = F, xlim = c(min(x.ticks), max(x.ticks)), xaxs = "i", ylim = c(0,1700), yaxs = "i")
fill.under.lines(X = snv[ ,Date], Y = snv[ ,New], YAxisMin = 0, Color = adjustcolor("black",.2))
abline(v = change.date, col = "black", xpd = F)
lines(x = snv[ ,Date], y = snv[ ,New], lwd = .7, col = adjustcolor("black",.5))
points(x = snv[ ,Date], y = snv[ ,New], pch = 16, col = adjustcolor(snv[ ,color.year.step],.7))
box()
axis(side = 1, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03, xpd = NA)
# axis(side = 1, at = x.lab.locs, labels = x.labs, lwd = 0, line = -.55, xpd = NA)
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 2, at = y.ticks, labels = T, lwd = 0, las = 2, line = -.3)
mtext(text = "New SNVs (count)", side = 2, line = 3.2, at = 800)

