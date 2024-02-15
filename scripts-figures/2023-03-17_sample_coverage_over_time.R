# RRR
# how well do the bins represent the samples over time?
# how well do the bins represent the taxonomic diversity of the lake?

# set-up ----

library(data.table)
library(lubridate)
source(file = "pop/scripts-generic_plots/plot_sample_dates_fancy.R")

yearly <- readRDS("pop/data/environmental_data/Robin-Refined/seasons/7-seasons_by_year.rds")
samply <- readRDS("pop/data/environmental_data/Robin-Refined/seasons/7-seasons_by_sample.rds")
bigtyme <- readRDS("data/2022-06-27_binning_groups/2022-06-27_bigtyme.rds")
color.key <- readRDS("pop/data/limony-IRD/2021-12-20_nmds_objects/Color_Keys.rds")
coverm <- readRDS("data/2023-03-15_coverM_on_drep96/processed/reads_mapped.rds")
all.bins <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
rel.abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")

options(scipen = 9999)

# ---- parse ----

yearly <- yearly[yearly$Year > 1999 & yearly$Year < 2020, ]

bigtyme <- merge(x = bigtyme, y = coverm$ANI_93, by.x = "sample.names", by.y = "sample", all = T)

colnames(rel.abund)[3] <- "abund.perc"
rel.abund <- rel.abund[Genome != "unmapped", .(Tot.Abund = sum(abund.perc)), by = .(Genome)]
tot.tot <- sum(rel.abund[ ,Tot.Abund])
rel.abund[ ,Ave.Abund := Tot.Abund / tot.tot * 100]
all.bins <- as.data.table(all.bins)
all.bins <- merge(x = all.bins, y = rel.abund[ ,.(Genome, Ave.Abund.perc = Ave.Abund)], by.x = "bin.full.name", by.y = "Genome", all.x = FALSE, all.y = TRUE)
all.bins[phylum == "p__" | is.na(phylum), phylum := "Unclassified"]
all.bins <- all.bins[order(Ave.Abund.perc)]

bin.count <- all.bins[ , .(.N, Tot.abund = sum(Ave.Abund.perc)),by = phylum]
bin.count <- bin.count[order(Tot.abund)]

all.bins <- merge(x = all.bins, y = bin.count, by = "phylum")
all.bins[Tot.abund < .2, phylum := "Other"]

bin.count <- all.bins[ , .(.N, Tot.abund = sum(Ave.Abund.perc)),by = phylum]
bin.count <- bin.count[order(Tot.abund)]

bin.stats <- all.bins[ ,.(bin.full.name, contamination, completeness)]

all.bins <- dcast(data = all.bins, formula = bin.full.name ~ phylum, value.var = "Ave.Abund.perc")
save.rows <- all.bins[ ,bin.full.name]
all.bins <- as.matrix(all.bins[ ,-1])
row.names(all.bins) <- save.rows
all.bins[is.na(all.bins)] <- 0
all.bins <- all.bins[order(rowSums(all.bins), decreasing = T),order(colSums(all.bins), decreasing = F)]
colnames(all.bins) <- sub(pattern = "p__", replacement = "", colnames(all.bins))
all.bins[1:5,1:5]



# plot ----

pdf(file = "figures/2023-03-17_fig_1_ideas/fig_1.pdf", width = 7.24, height = 2.24)

# season shading 

par(mar = c(1.75,2.6,.1,5.5), fig = c(.4,1,0,1))
y.ax.labels <- set.up.plot(dates.vector = bigtyme$sample.dates, year.lines = F, month.lines = F, 
                           x.ax.line = 0, x.ax.lab.line = -1.15, x.ax.single.letters = T,
                           y.ax.line = 0, y.ax.lab.line = -0.6, y.ax.btwn.lines = F,
                           x.ax.cex = .7, y.ax.cex = .7)
mtext(text = "Sample Year", side = 2, line = 1.9, cex = .9)
mtext(text = "Sample Month", side = 1, line = .75, cex = .9)

shade.by.season(season.starts = yearly$Spring.Ice.On, season.ends = yearly$Spring.Ice.Off, season.color = adjustcolor("snow2",.7))
shade.by.season(season.starts = yearly$Fall.Ice.On, season.ends = yearly$Fall.Ice.Off, season.color =  adjustcolor("snow2",.7))

shade.by.season(season.starts = yearly$Spring.Ice.Off, season.ends = yearly$Clear.Start, season.color = adjustcolor("tan4",.7))

shade.by.season(season.starts = yearly$Clear.Start, season.ends = yearly$Clear.End, season.color = adjustcolor("cornflowerblue",.7))

shade.by.season(season.starts = yearly$Clear.End, season.ends = yearly$Date.Samples.First.Anoxic, season.color =  adjustcolor("chartreuse4",.7))

shade.by.season(season.starts = yearly$Date.Samples.First.Anoxic, season.ends = yearly$Date.Temp.Diff.Below.3, season.color = adjustcolor("purple",.3))

i.iceless.fall <- which(is.na(yearly$Fall.Ice.On))
yearly$Fall.Ice.On[i.iceless.fall] <- parse_date_time(paste(yearly$Year[i.iceless.fall],"12 31"), orders = "ymd", tz = "Etc/GMT-5")
shade.by.season(season.starts = yearly$Date.Temp.Diff.Below.3, season.ends = yearly$Fall.Ice.On, season.color = adjustcolor("hotpink2",.5))
i.iceless.winter <- which(yearly$Spring.Ice.On != parse_date_time(paste(yearly$Year,1,1), orders = "ymd"))
shade.by.season(season.starts = parse_date_time(paste(yearly$Year[i.iceless.winter],1,1), orders = "ymd", tz = "Etc/GMT-5"), 
                season.ends = yearly$Spring.Ice.On[i.iceless.winter], season.color = adjustcolor("hotpink2",.5))
box()

# dates scaled by bin coverage 

bigtyme$cex <- get.scaled.point.cex(min.cex = .2, max.cex = 1, my.values = bigtyme$perc.mapped)
summary(bigtyme$perc.mapped)
legend.vals <- c(10,20,30,40,50,60)[6:1]
legend.cex <- get.scaled.point.cex(min.cex = .2, max.cex = 1, my.values = legend.vals)

add.dates(date.vector = bigtyme$sample.dates, year.labs = y.ax.labels, 
          point.type = 21, point.col = adjustcolor("black", alpha.f = .2), line.col = adjustcolor("black",.7), 
          point.cex = bigtyme$cex)

# top legend
# add.legend.point(date.vector = c("8-15-2020","8-15-2021"), point.vertical.scoot = .25,
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = .3), line.col = adjustcolor("black",.9),
#                  point.cex = legend.cex[1], description = "", text.adj = .5, text.cex = .8)
# add.legend.point(date.vector = c("9-1-2020","9-1-2021"), point.vertical.scoot = .25,
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = .3), line.col = adjustcolor("black",.9),
#                  point.cex = legend.cex[2], description = legend.vals[2], text.adj = .5, text.cex = .8)
# add.legend.point(date.vector = c("9-15-2020","9-15-2021"), point.vertical.scoot = .25,
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = .3), line.col = adjustcolor("black",.9),
#                  point.cex = legend.cex[3], description = "", text.adj = .5, text.cex = .8)
# add.legend.point(date.vector = c("10-1-2020","10-1-2021"), point.vertical.scoot = .25,
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = .3), line.col = adjustcolor("black",.9),
#                  point.cex = legend.cex[4], description = legend.vals[4], text.adj = .5, text.cex = .8)
# add.legend.point(date.vector = c("10-15-2020","10-15-2021"), point.vertical.scoot = .25,
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = .3), line.col = adjustcolor("black",.9),
#                  point.cex = legend.cex[5], description = "", text.adj = .5, text.cex = .8)
# add.legend.point(date.vector = c("11-1-2020","11-1-2021"), point.vertical.scoot = .25,
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = .3), line.col = adjustcolor("black",.9),
#                  point.cex = legend.cex[6], description = legend.vals[6], text.adj = .5, text.cex = .8)
# add.legend.point(date.vector = c("8-1-2020","8-1-2020"), text.vertical.scoot = .5,
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = 0), line.col = adjustcolor("black",0),
#                  description = "Reads\nMapped (%)", text.adj = 1)
# rect(xleft = parse_date_time("5-15-2000","mdy"), xright = parse_date_time("11-10-2000","mdy"),
#      ybottom = 2019.75, ytop = 2021.5, xpd = NA)

# side legend
add.side.legend.points(year = seq(from = 2015.5, to = 2010, along.with = legend.vals), month.day = "1-30",
                       point.type = 21, point.col = adjustcolor("black", alpha.f = .2), line.col = adjustcolor("black",.7),
                       point.cex = legend.cex)
add.side.legend.text(year = seq(from = 2015.5, to = 2010, along.with = legend.vals), month.day = "2-15", 
                     description = legend.vals, text.adj = 0, text.cex = .7, text.col = c("black",adjustcolor("black",alpha.f = 0)))
add.side.legend.text(year = 2018.9, month.day = "1-20", description = "Reads", text.cex = .9)
add.side.legend.text(year = 2017.3, month.day = "1-20", description = "Mapped (%)", text.cex = .9)

add.side.legend.text(year = 2007.25, month.day = "1-20", description = "Season", text.cex = .9)
add.side.legend.points(year = seq(from = 2005.5, to = 1997, along.with = legend.vals), month.day = "1-30",
                       point.type = 22, point.cex = 1.5, line.col = "black",
                       point.col = c(adjustcolor("snow2",.7), adjustcolor("tan4",.7), adjustcolor("cornflowerblue",.7), adjustcolor("chartreuse4",.7), adjustcolor("purple",.3), adjustcolor("hotpink2",.5)))
add.side.legend.text(year = seq(from = 2005.5, to = 1997, along.with = legend.vals), month.day = "2-15", 
                     description = c("Ice-On","Spring","Clearwater","Early Summer","Late Summer","Fall"), text.adj = 0, text.cex = .7)



# # dates shaded by bin coverage - looks bad
# 
# bigtyme$alpha <- bigtyme$perc.mapped / 100 
# bigtyme$color <- "black"
# for(r in 1:nrow(bigtyme)){
#   bigtyme$color[r] <- adjustcolor("black", alpha.f = bigtyme$alpha[r])
# }
# 
# legend.cex <- pi * (c(20,50,80) / 80)^2
# 
# add.dates(date.vector = bigtyme$sample.dates, year.labs = y.ax.labels, 
#           point.type = 21, point.col = bigtyme$color, line.col = adjustcolor("black",1), 
#           point.cex = 1)
# add.legend.point(date.vector = c("3-1-2020","3-1-2021"), point.vertical.scoot = .25, 
#                  point.type = 21, point.col = bigtyme$color, line.col = adjustcolor("black",.9), 
#                  point.cex = legend.cex[1], description = "20%", text.adj = .5)
# add.legend.point(date.vector = c("4-1-2020","4-1-2021"), point.vertical.scoot = .25, 
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = .3), line.col = adjustcolor("black",.9), 
#                  point.cex = legend.cex[2], description = "50%", text.adj = .5)
# add.legend.point(date.vector = c("5-1-2020","5-1-2021"), point.vertical.scoot = .25, 
#                  point.type = 21, point.col = adjustcolor("black", alpha.f = .3), line.col = adjustcolor("black",.9), 
#                  point.cex = legend.cex[3], description = "80%", text.adj = .5)

# ---- make rank abund plot 

par(mar =c(2.1,5,.2,2), fig= c(0,.4,0,1), new = T)

bar.spots <- barplot(height = all.bins, beside = F, horiz = T, ann = F, axes = F, col = adjustcolor("black", alpha.f = .3), names.arg = rep("",ncol(all.bins)), space = .5, xaxs = "i", yaxs = "i")
text(x = -2, y = bar.spots, labels = colnames(all.bins), xpd = NA, cex = .7, adj = 1)
axis(side = 1, labels = F, line = .25, tck = -.04, at = c(0,10,20,30,40), xpd = NA)
axis(side = 1, labels = T, lwd = 0, line = -.75, cex.axis = .7, at = c(0,10,20,30,40), xpd = NA)
text(x = bin.count$Tot.abund + 2, y = bar.spots, labels = bin.count$N, cex = .7, adj = 0, xpd = NA)
mtext(text = "Abundance (%)", side = 1, line = 1.1, cex = .9)

rect(xleft = 14, xright = 42, ybottom = 15, ytop = 1, xpd = NA)

# ---- make completeness/quality plot 

par(mar =c(0,0,0,0), fig = c(.265,.395,.265,.665), new = T)
plot(x = c(0,5), y = c(50,100), type = "n", axes = F, ann = F)
# box(which = "figure")
boxplot(x = bin.stats[ ,.(completeness)], add = T, ann = F, axes = F, at = 1, lty = 1, range = 0, col = adjustcolor("blue",.5))
axis(side = 2, at = c(50,75,100), labels = F, tck = -.025, line = -.6, col = "blue")
axis(side = 2, at = c(50,75,100), cex.axis = .7, labels = T, lwd = 0, line = -1.3, las = 2, col.axis = "blue")
mtext(text = "Completeness", side = 2, line = .5, cex = .7, col = "blue")

par(mar =c(0,0,0,0), fig = c(.265,.395,.265,.665), new = T)
plot(x = c(0,5), y = c(0,10), type = "n", axes = F, ann = F)
# box(which = "figure")
boxplot(x = bin.stats[ ,.(contamination)], add = T, ann = F, axes = F, at = 2, lty = 1, range = 0, xpd = NA, col = adjustcolor("red4",.5))
axis(side = 4, at = c(0,5,10), labels = F, tck = -.025, line = -2.3, col = "red4")
axis(side = 4, at = c(0,5,10), cex.axis = .7, labels = T, lwd = 0, line = -3, las = 2, col.axis = "red4")
mtext(text = "Contamination", side = 4, line = -1.9, cex = .7, col = "red4")


# 

dev.off()

# ----
