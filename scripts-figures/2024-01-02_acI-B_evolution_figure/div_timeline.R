# RRR
# data for the acI-B figure: ME2011-09-21_3300043464_group3_bin69
# Fig 5b - diversity over time

library(data.table)
library(lubridate)

div <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/abund_and_nucl_div-ME2011-09-21_3300043464_group3_bin69.csv", colClasses = c("date" = "character"))

index.change <- which(div$date == "2012-08-03")

div[ ,date := parse_date_time(date, "ymd")]

change.date <- div$date[index.change]

fill.under.lines <- function(X, Y, YAxisMin, Color, xpd = F){
  poly.x <- c(min(X), X, max(X))
  poly.y <- c(YAxisMin, Y, YAxisMin )
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

x.ticks <- parse_date_time(paste(seq(2000,2020,2), 1, 1), "ymd")
# x.labs <- year(x.ticks)
# x.lab.locs <- x.ticks
# x.lab.locs[x.labs == 2020] <- x.lab.locs[x.labs == 2020] %m+% -months(0)

summary(div$nucl_diversity)
y.ticks <- seq(.008,.018,.003)

plot(x = div[ ,date], y = div[ ,nucl_diversity], type = "n", ann = F, axes = F, xlim = c(min(x.ticks), max(x.ticks)), xaxs = "i")
fill.under.lines(X = div[ ,date], Y = div[ ,nucl_diversity], YAxisMin = 0, Color = adjustcolor("grey40",.2))
abline(v = change.date, col = "black", xpd = F)
lines(x = div[ ,date], y = div[ ,nucl_diversity], lwd = .7, col = adjustcolor("black",.5))
points(x = div[ ,date], y = div[ ,nucl_diversity], pch = 16, col = adjustcolor(div[ ,color.year.step], .7))
box()
axis(side = 1, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03, xpd = NA)
# axis(side = 1, at = x.lab.locs, labels = x.labs, lwd = 0, line = -.55, xpd = NA)
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 2, at = y.ticks, labels = T, lwd = 0, las = 2, line = -.3)
mtext(text = expression(paste("Diversity (",pi,")")), side = 2, line = 3)

