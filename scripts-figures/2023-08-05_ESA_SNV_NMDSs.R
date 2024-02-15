# RRR
# Make NMDS plots with base R, force them to be square etc
# make sure lines are actually going in order!
# Use Euclidean distance because total SNVs does not have a max carrying capacity, and increases after 2012

library(data.table)
library(lubridate)
library(viridisLite)

my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-B_nmds_euc.rds")
my.genome <- "acI-B"

my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-A_nmds_euc.rds")
my.genome <- "acI-A"

my.nmds <- readRDS(file = "figures/2023-07-21_NMDS_plots_of_SNVs/acI-C_nmds_euc.rds")
my.genome <- "acI-C"

# ---- set up for plot ----

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
pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-B_long-term.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
my.lims <- c(-max(c(my.nmds$x, my.nmds$y))-5, max(c(my.nmds$x, my.nmds$y))+5)
plot(x = my.nmds$x, y = my.nmds$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
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
pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-B_long-term-1.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
my.lims <- c(-max(c(my.nmds$x, my.nmds$y))-5, max(c(my.nmds$x, my.nmds$y))+5)
plot(x = my.nmds$x, y = my.nmds$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
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

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-B_long-term-2.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
my.lims <- c(-max(c(my.nmds$x, my.nmds$y))-5, max(c(my.nmds$x, my.nmds$y))+5)
plot(x = my.nmds$x, y = my.nmds$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
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

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-B_long-term-3.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
my.lims <- c(-max(c(my.nmds$x, my.nmds$y))-5, max(c(my.nmds$x, my.nmds$y))+5)
plot(x = my.nmds$x, y = my.nmds$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
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

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-B_long-term-4.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
my.lims <- c(-max(c(my.nmds$x, my.nmds$y))-5, max(c(my.nmds$x, my.nmds$y))+5)
plot(x = my.nmds$x, y = my.nmds$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
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

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-B_long-term-5.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
my.lims <- c(-max(c(my.nmds$x, my.nmds$y))-5, max(c(my.nmds$x, my.nmds$y))+5)
plot(x = my.nmds$x, y = my.nmds$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
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

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-B_long-term-6.pdf", width = 4.2, height = 3)
par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
my.lims <- c(-max(c(my.nmds$x, my.nmds$y))-5, max(c(my.nmds$x, my.nmds$y))+5)
plot(x = my.nmds$x, y = my.nmds$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
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

# ----

p.drift <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = year), alpha = .7, size = 3)+
  scale_color_viridis_c(option = "D")+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  guides(color = guide_legend(title = "Year"))

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-C_drift.pdf", width = 5, height = 4)
p.drift
dev.off()

# ----

# ICE
my.colors <- c("snow3","wheat","wheat","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(4,1.5,1.5,1.5,1.5,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.ice <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Ice-On")

# SPRING
my.colors <- c("wheat","tan4","wheat","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,4,1.5,1.5,1.5,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.spring <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),        
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Spring")

# CLEAR
my.colors <- c("wheat","wheat","cornflowerblue","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,1.5,4,1.5,1.5,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.clear <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Clearwater")

# EARLY
my.colors <- c("wheat","wheat","wheat","chartreuse4","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,1.5,1.5,4,1.5,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.early <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Early Summer")

# LATE
my.colors <- c("wheat","wheat","wheat","wheat","purple","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,1.5,1.5,1.5,4,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.late <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Late Summer")

# FALL
my.colors <- c("wheat","wheat","wheat","wheat","wheat","hotpink2")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,1.5,1.5,1.5,1.5,4)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.fall <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Fall")

p.seasons.separate <- p.spring + p.clear + p.early + p.late + p.fall + p.ice + plot_layout(guides = "collect")

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-C_seasons.pdf", width = 10, height = 7.25)
p.seasons.separate
dev.off()
# ----

# acI-A NMDS seasons plots ----

library(ggplot2)
library(patchwork)

my.nmds[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]

# ----

p.drift <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_point(aes(color = year), alpha = .7, size = 3)+
  scale_color_viridis_c(option = "D")+
  ylab(label = "NMDS Axis 2")+
  xlab(label = "NMDS Axis 1")+
  guides(color = guide_legend(title = "Year"))

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-A_drift.pdf", width = 5, height = 4)
p.drift
dev.off()

# ----

# ICE
my.colors <- c("snow3","wheat","wheat","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(4,1.5,1.5,1.5,1.5,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.ice <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Ice-On")

# SPRING
my.colors <- c("wheat","tan4","wheat","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,4,1.5,1.5,1.5,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.spring <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),        
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Spring")

# CLEAR
my.colors <- c("wheat","wheat","cornflowerblue","wheat","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,1.5,4,1.5,1.5,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.clear <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Clearwater")

# EARLY
my.colors <- c("wheat","wheat","wheat","chartreuse4","wheat","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,1.5,1.5,4,1.5,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.early <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Early Summer")

# LATE
my.colors <- c("wheat","wheat","wheat","wheat","purple","wheat")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,1.5,1.5,1.5,4,1.5)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.late <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Late Summer")

# FALL
my.colors <- c("wheat","wheat","wheat","wheat","wheat","hotpink2")
names(my.colors) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
my.size <- c(1.5,1.5,1.5,1.5,1.5,4)
names(my.size) <- c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")
p.fall <- ggplot(data = my.nmds, aes(x = x, y = y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 25))+
  geom_point(aes(color = season, size = season), alpha = .7, show.legend = F)+
  geom_path(alpha = 0)+
  scale_color_manual(values = my.colors)+
  scale_size_manual(values = my.size)+
  labs(title = "Fall")

p.seasons.separate <- p.spring + p.clear + p.early + p.late + p.fall + p.ice + plot_layout(guides = "collect")

pdf(file = "figures/2023-08-05_ESA_plots/NMDS_acI-A_seasons.pdf", width = 10, height = 7.25)
p.seasons.separate
dev.off()
# ----
