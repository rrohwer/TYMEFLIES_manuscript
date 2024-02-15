# RRR
# data for the acI-B figure: ME2011-09-21_3300043464_group3_bin69
# Fig 5a - abundance over time

library(data.table)
library(lubridate)

abund <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/abund_and_nucl_div-ME2011-09-21_3300043464_group3_bin69.csv", colClasses = c("date" = "character"))

index.change <- which(abund$date == "2012-08-03")

abund[ ,date := parse_date_time(date, "ymd")]

change.date <- abund$date[index.change]

fill.under.lines <- function(X, Y, YAxisMin, Color, xpd = F){
  poly.x <- c(min(X), X, max(X))
  poly.y <- c(YAxisMin, Y, YAxisMin )
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

x.ticks <- parse_date_time(paste(seq(2000,2020,2), 1, 1), "ymd")
# x.labs <- year(x.ticks)
# x.lab.locs <- x.ticks
# x.lab.locs[x.labs == 2020] <- x.lab.locs[x.labs == 2020] %m+% -months(0)

summary(abund$adj.abund.perc)
y.ticks <- seq(0,.7,.2)

plot(x = abund[ ,date], y = abund[ ,adj.abund.perc], type = "n", ann = F, axes = F, xlim = c(min(x.ticks), max(x.ticks)), xaxs = "i")
fill.under.lines(X = abund[ ,date], Y = abund[ ,adj.abund.perc], YAxisMin = 0, Color = adjustcolor("grey40",.2))
abline(v = change.date, col = "black", xpd = F)
lines(x = abund[ ,date], y = abund[ ,adj.abund.perc], lwd = .7, col = adjustcolor("black",.5))
points(x = abund[ ,date], y = abund[ ,adj.abund.perc], pch = 16, col = adjustcolor(abund[ ,color.year.step],.7))
box()
axis(side = 1, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03, xpd = NA)
# axis(side = 1, at = x.lab.locs, labels = x.labs, lwd = 0, line = -.55, xpd = NA)
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 2, at = y.ticks, labels = T, lwd = 0, las = 2, line = -.3)
mtext(text = "Abundance (%)", side = 2, line = 3.2)
text(x = change.date %m+% -days(90), y = .68, labels = "2012-08-03", adj = 1)



