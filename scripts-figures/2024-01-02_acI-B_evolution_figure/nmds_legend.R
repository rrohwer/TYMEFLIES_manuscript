# RRR
# draw the legend separately. adjust nmds dims to be exactly square, legend has more give


# # this was the vertical one:
# 
# y.locs <- seq(from = 8, to = 4, along.with = 1:6)
# y.locs[1] <- y.locs[1] - (y.locs[2]-y.locs[1]) / 4
# y.locs[6] <- y.locs[6] + (y.locs[2]-y.locs[1]) / 4
# title.loc <- y.locs[1] - (y.locs[2]-y.locs[1]) 
# 
# x.points <- seq(from = 0, to = 0, along.with = 1:6)
# x.text <- seq(from = 2, to = 2, along.with = 1:6)
# 
# plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
# 
# points(y = y.locs, x = x.points, pch = 16, cex = 2, xpd = NA,
#        col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)))
# 
# text(y = y.locs, x = x.text, labels = c("2000 -\n2010", "2011","2012","2013","2014","2015 -\n2019"), adj = 0, xpd = NA)
# 
# text(x = 0, y = title.loc, labels = "Year:", adj = 0, xpd = NA)

# but make it horizontal

x.points <- seq(from = 2, to = 8, along.with = 1:6)
title.loc <- x.points[1] - (x.points[2]-x.points[1]) /1.25

x.text <- x.points + .25

y.locs <- seq(from = 5, to = 5, along.with = 1:6)

plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)

points(x = x.points, y = y.locs, pch = 16, cex = 2, xpd = NA,
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)))

text(x = x.text, y = y.locs, labels = c("2000 -\n2010", "2011","2012","2013","2014","2015 -\n2019"), adj = 0, xpd = NA)

text(x = title.loc, y = y.locs[1], labels = "Year:", adj = 0, xpd = NA)
