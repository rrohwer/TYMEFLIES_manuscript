# RRR
# Make NMDS plots with base R, force them to be square etc
# make sure lines are actually going in order!
# Use Euclidean distance because total SNVs does not have a max carrying capacity, and increases after 2012

library(data.table)
library(lubridate)
library(viridisLite)

sample.key <- fread(file = "data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv")

my.nmds <- readRDS(file = "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-21_3300043464_group3_bin69_all_SNV_euclidean_distance_NMDS_object.rds")
my.genome <- "acI-B"

# my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-A_nmds_euc.rds")
# my.genome <- "acI-A"

my.nmds <- readRDS(file = "data/2023-11-01_multidimensional_SNV_analysis/ME2016-07-20_3300033996_group7_bin32_all_SNV_euclidean_distance_NMDS_object.rds")
my.genome <- "acI-C"

# ---- set up for plot ----

my.nmds <- data.table("sample" = rownames(my.nmds$points),"x" = my.nmds$points[ ,1], "y" = my.nmds$points[ ,2])
my.nmds <- merge(x = my.nmds, y = sample.key, by = "sample")
my.nmds[ ,`:=`(season = factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")),
               invasion = factor(invasion, levels = c("none","spiny","zebra")),
               date = parse_date_time(x = date, orders = "ymd"))]
my.nmds <- my.nmds[order(date)]


color.key <- data.table("season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"),
                        "color.season" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))
my.nmds <- merge(x = my.nmds, y = color.key, by = "season")

color.key <- data.table("year" = 2000:2019, 
                        "color.year.gradual" = viridis(n = 20, option = "B"))
my.nmds <- merge(x = my.nmds, y = color.key, by = "year")               

color.key <- data.table("year" = 2000:2019, 
                        "color.year.step" = c(rep("#b89996", 11), "purple","red","tan2","dodgerblue4", rep("#7592b2", 5)))
my.nmds <- merge(x = my.nmds, y = color.key, by = "year")               

background.index <- which(my.nmds$year != 2011 & my.nmds$year != 2012 & my.nmds$year != 2013 & my.nmds$year != 2014)
forground.index  <- which(my.nmds$year == 2011 | my.nmds$year == 2013 | my.nmds$year == 2014)
index.2012 <- which(my.nmds$year == 2012)

# ---- make plots ----

# acI-B long-term step change ----
pdf(file = "figures/2023-11-07_scripps_plots/NMDS_acI-B_long-term.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
x.lim <- c(-35,55) # range of 90
y.lim <- c(-44,46) # range of 90
plot(x = my.nmds$x, y = my.nmds$y, xlim = x.lim, ylim = y.lim, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = my.nmds$x[background.index], y = my.nmds$y[background.index], pch = 16, col = adjustcolor(my.nmds$color.year.step[background.index], alpha.f = .7))
points(x = my.nmds$x[forground.index], y = my.nmds$y[forground.index], pch = 16, col = adjustcolor(my.nmds$color.year.step[forground.index], alpha.f = .7))
points(x = my.nmds$x[index.2012], y = my.nmds$y[index.2012], pch = 16, col = adjustcolor(my.nmds$color.year.step[index.2012], alpha.f = 1))
par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 2, to = 2, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)

dev.off()
# ----

# acI-B animated ----
pdf(file = "figures/2023-11-07_scripps_plots/NMDS_acI-B_long-term-1.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(smar = c(1.5,1.5,1.5,1.5))
x.lim <- c(-35,55) # range of 90
y.lim <- c(-44,46) # range of 90
plot(x = my.nmds$x, y = my.nmds$y, xlim = x.lim, ylim = y.lim, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = my.nmds$x[my.nmds$year < 2011], y = my.nmds$y[my.nmds$year < 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year < 2011], alpha.f = .7))
par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 2, to = 2, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)
dev.off()

pdf(file = "figures/2023-11-07_scripps_plots/NMDS_acI-B_long-term-2.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = my.nmds$x, y = my.nmds$y, xlim = x.lim, ylim = y.lim, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = my.nmds$x[my.nmds$year < 2011], y = my.nmds$y[my.nmds$year < 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year < 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2011], y = my.nmds$y[my.nmds$year == 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2011], alpha.f = .7))
par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 2, to = 2, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)
dev.off()

pdf(file = "figures/2023-11-07_scripps_plots/NMDS_acI-B_long-term-3.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = my.nmds$x, y = my.nmds$y, xlim = x.lim, ylim = y.lim, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = my.nmds$x[my.nmds$year < 2011], y = my.nmds$y[my.nmds$year < 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year < 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2011], y = my.nmds$y[my.nmds$year == 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2012], y = my.nmds$y[my.nmds$year == 2012], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2012], alpha.f = 1))
par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 2, to = 2, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)
dev.off()

pdf(file = "figures/2023-11-07_scripps_plots/NMDS_acI-B_long-term-4.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = my.nmds$x, y = my.nmds$y, xlim = x.lim, ylim = y.lim, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = my.nmds$x[my.nmds$year < 2011], y = my.nmds$y[my.nmds$year < 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year < 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2011], y = my.nmds$y[my.nmds$year == 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2012], y = my.nmds$y[my.nmds$year == 2012], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2012], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2013], y = my.nmds$y[my.nmds$year == 2013], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2013], alpha.f = .7))
par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 2, to = 2, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)
dev.off()

pdf(file = "figures/2023-11-07_scripps_plots/NMDS_acI-B_long-term-5.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = my.nmds$x, y = my.nmds$y, xlim = x.lim, ylim = y.lim, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = my.nmds$x[my.nmds$year < 2011], y = my.nmds$y[my.nmds$year < 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year < 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2011], y = my.nmds$y[my.nmds$year == 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2012], y = my.nmds$y[my.nmds$year == 2012], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2012], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2013], y = my.nmds$y[my.nmds$year == 2013], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2013], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2014], y = my.nmds$y[my.nmds$year == 2014], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2014], alpha.f = .7))
par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 2, to = 2, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)
dev.off()

pdf(file = "figures/2023-11-07_scripps_plots/NMDS_acI-B_long-term-6.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = my.nmds$x, y = my.nmds$y, xlim = x.lim, ylim = y.lim, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = my.nmds$x[my.nmds$year < 2011], y = my.nmds$y[my.nmds$year < 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year < 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2011], y = my.nmds$y[my.nmds$year == 2011], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2011], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2012], y = my.nmds$y[my.nmds$year == 2012], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2012], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2013], y = my.nmds$y[my.nmds$year == 2013], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2013], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year == 2014], y = my.nmds$y[my.nmds$year == 2014], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year == 2014], alpha.f = .4))
points(x = my.nmds$x[my.nmds$year > 2014], y = my.nmds$y[my.nmds$year > 2014], pch = 16, col = adjustcolor(my.nmds$color.year.step[my.nmds$year > 2014], alpha.f = .7))
par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 2, to = 2, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)
dev.off()

# ----

# acI-C NMDS seasons plots ----

library(ggplot2)
library(patchwork)

my.nmds[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
my.nmds[ ,year := as.numeric(year)]
my.nmds[ ,Year := year]

# ----

p.drift <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = Year), alpha = .7, size = 3)+
  scale_color_viridis_c(option = "viridis")+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")

pdf(file = "figures/2023-11-07_scripps_plots/NMDS_acI-C_drift.pdf", width = 5, height = 4)
p.drift
dev.off()

# ----

