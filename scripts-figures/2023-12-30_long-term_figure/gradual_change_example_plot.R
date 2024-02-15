# RRR
# fig 3a

library(data.table)
library(lubridate)

dist.first <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/snv_dist_table-gradual-ME2005-06-22_3300042363_group2_bin84-.csv", colClasses = c("date.1" = "character","date.2" = "character"))

dist.first[ ,`:=`(date.1 = parse_date_time(date.1, "ymd"), date.2 = parse_date_time(date.2, "ymd"))]


x.ticks <- parse_date_time(paste(seq(2000,2020,2), 1, 1), "ymd")
x.labs <- year(x.ticks)
x.lab.locs <- x.ticks
x.lab.locs[x.labs == 2020] <- x.lab.locs[x.labs == 2020] %m+% -months(4)

y.ticks <- seq(25,60,15)

plot(x = dist.first[ ,date.2], y = dist.first[ ,dist], type = "n", ann = F, axes = F, xlim = c(min(x.ticks),max(x.ticks)), xaxs = "i")
lines(x = dist.first[ ,date.2], y = dist.first[ ,dist], lwd = .7, col = adjustcolor("black",.5))
points(x = dist.first[ ,date.2], y = dist.first[ ,dist], pch = 16, col = adjustcolor("black",.7))
box()
axis(side = 1, at = x.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 1, at = x.lab.locs, labels = x.labs, lwd = 0, line = -.55)
axis(side = 2, at = y.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 2, at = y.ticks, labels = T, lwd = 0, las = 2, line = -.3)
mtext(text = "Distance", side = 2, line = 2)
mtext(text = "Gradual Change", side = 3, at = parse_date_time(x = paste(2000,10,1), orders = "ymd"), line = -2, adj = 0)
