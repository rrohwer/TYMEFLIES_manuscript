# RRR
# fig 4e
# yahara discharge
# include 2023 because the year is almost over 

library(data.table)

discharge <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/yahara_discharge.csv")

x.ax.ticks <- yday(parse_date_time(paste(1:12,"-1"), orders = "md"))[seq.int(1,12,2)]
x.ax.lab.locs <- x.ax.ticks # yday(parse_date_time(paste(1:12,"-15"), orders = "md"))
x.ax.labs <- month(parse_date_time(paste(1:12,"-1"), orders = "md"), label = TRUE, abbr = TRUE)[seq.int(1,12,2)]

# change units so that it takes less space on the y-axis

discharge[ ,discharge.m3.s := discharge.ft3.s * 0.0283168]

y.lim <- c(min(discharge$discharge.m3.s), max(discharge$discharge.m3.s))
y.ax.ticks <- seq.int(0,36,9)


plot(x = c(0,365), y = y.lim, type = "n", axes = F, ann = F)
box()
axis(side = 1, at = x.ax.ticks, labels = F, lwd = 0, lwd.ticks = 1, tck = -.03)
axis(side = 1, at = x.ax.lab.locs[seq.int(1,12,2)], labels = x.ax.labs[seq.int(1,12,2)], lwd = 0, las = 1, line = -.6)
axis(side = 1, at = x.ax.lab.locs[seq.int(2,12,2)], labels = x.ax.labs[seq.int(2,12,2)], lwd = 0, las = 1, line = -.6)
axis(side = 2, lwd = 0, at = y.ax.ticks, lwd.ticks = 1, labels = F, tck = -.03)
axis(side = 2, lwd = 0, at = y.ax.ticks, labels = T, las = 2, line = -.5)
mtext(text = expression(Discharge~(m^3/s)), side = 2, line = 1.8)
# for those fucking superscripts! https://www.dataanalytics.org.uk/axis-labels-in-r-plots-using-expression/
for(yr in unique(discharge$year)){
  lines(x = discharge[year == yr, yday], y = discharge[year == yr, discharge.m3.s], col = adjustcolor("grey40",.4), lwd = .7)
}
lines(x = discharge[year == 2012, yday], y = discharge[year == 2012, discharge.m3.s], col = "red2", lwd = 1)
