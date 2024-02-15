# RRR

# make figure 4b
# epilimnion mean temperature in 2012 

library(data.table)

temp <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/epilimnion_mean_temperatures.csv")

# ---- make plot ----

x.ax.ticks <- yday(parse_date_time(paste(1:12,"-1"), orders = "md"))[seq.int(1,12,2)]
x.ax.lab.locs <- x.ax.ticks # yday(parse_date_time(paste(1:12,"-15"), orders = "md"))
x.ax.labs <- month(parse_date_time(paste(1:12,"-1"), orders = "md"), label = TRUE, abbr = TRUE)[seq.int(1,12,2)]

y.ax.ticks <- seq.int(0,30,10)

plot(x = c(0,365), y = c(min(temp[ , Temp.C], na.rm = T), max(temp[ , Temp.C], na.rm = T)), type = "n", axes = F, ann = F)
box()
axis(side = 1, at = x.ax.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 1, at = x.ax.lab.locs[seq.int(1,12,2)], labels = x.ax.labs[seq.int(1,12,2)], lwd = 0, las = 1, line = -.6)
axis(side = 1, at = x.ax.lab.locs[seq.int(2,12,2)], labels = x.ax.labs[seq.int(2,12,2)], lwd = 0, las = 1, line = -.6)
axis(side = 2, lwd = 0, at = y.ax.ticks, lwd.ticks = 1, labels = F, tck = -.03)
axis(side = 2, lwd = 0, at = y.ax.ticks, labels = T, las = 2, line = -.5)
mtext(text = "Temperature (°C)", side = 2, line = 2)
# mtext(text = "Month", side = 1, line = 2)

for (yr in unique(temp[Year != 2012, Year])){
  points(x = temp[Year == yr, yDay], y = temp[Year == yr, Temp.C], col = adjustcolor("grey40",.2), pch = 19, cex = .1)
}
lines(x = temp[Year == 2012, yDay], y = temp[Year == 2012, Temp.C], col = adjustcolor("red2",1), lwd = 1)

segments(x0 = 10, x1 = 30, y0 = 25, y1 = 25, col = "red2", lwd = 2)
text(x = 40, y = 25, label = "2012", adj = 0)
# points(x = 20, y = 20, col = "red2", pch = 19, cex = 1)
# text(x = 40, y = 20, label = "2012", adj = 0)
