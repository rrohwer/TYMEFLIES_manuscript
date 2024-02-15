# RRR
# Make the abundance over time plot a single color

library(data.table)
library(lubridate)

coverm <- fread(input = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
genome.info <- fread(input = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")

my.genome <- "ME2011-09-21_3300043464_group3_bin69" # acI-B

breadth.cutoff <- .5

# ---- parse data ----

colnames(coverm)[3] <- "rel.abund"
coverm <- coverm[ ,1:3]
coverm <- coverm[Genome == my.genome]
coverm[ ,date := substr(Sample, start = 3, stop = 12)]
coverm[ ,date := parse_date_time(x = date, orders = "ymd")]
coverm[ ,Sample := NULL]
coverm <- coverm[ ,.(rel.abund = mean(rel.abund)), by = date] # should choose unfiltered for the one dup day with ww and pf

genome.info <- genome.info[genome == my.genome]
genome.info[ ,date := parse_date_time(x = date, orders = "ymd")]
genome.info[ ,above.breadth.cutoff := (breadth / breadth_expected) >= breadth.cutoff, ]
genome.info <- genome.info[ ,.(nucl_diversity = mean(nucl_diversity)), by = .(date, year, yday, season, invasion, above.breadth.cutoff)]

info <- merge(x = genome.info, y = coverm, by = "date")

color.key <- data.table("season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"),
                        "color" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))
info <- merge(x = info, y = color.key, by = "season")
info <- info[order(date)]

# ---- plot set-up base R ----

add.magic.boxes <- function(){
  box(which = "inner", col="red", lwd = 3)
  box(which = "outer", col="blue", lwd = 3)
  box(which = "plot", col="purple", lwd = 3)
  box(which = "figure", col="orange", lwd = 3)
  # inner margins (orange to purple) are mar/plt/mai
  # outer margins (red to blue) are oma/omd/omi
  par.settings <- par(no.readonly = TRUE)
  cat("change distance btwn ORANGE and PURPLE with mar  (mar = ", par.settings$mar, ")\n")
  cat("change distance btwn BLUE and RED with oma       (oma = ", par.settings$oma, ")\n")
}

fill.under.lines <- function(X, Y, YAxisMin, Color, xpd = F){
  poly.x <- c(min(X), X, max(X))
  poly.y <- c(YAxisMin, Y, YAxisMin )
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

x.ticks <- parse_date_time(x = c(2000,2012,2019), orders = "y")

# # acI-B abund ----
# pdf(file = "figures/2023-08-10_UCSB_plots/abund_acI-B_long-term.pdf",width = 6.5, height = 3)
# par(mar = c(1.75,3,1,0))
# plot(x = info$date, y = info$rel.abund, type = "n", ann = F, axes = F)
# fill.under.lines(X = info$date, Y = info$rel.abund, YAxisMin = 0, Color = adjustcolor(col = "#b89996", alpha.f = .5))
# lines(x = info$date, y = info$rel.abund, col = "#b89996")
# points(x = info$date, y = info$rel.abund, pch = 16, col = adjustcolor("#b89996", alpha.f = .7))
# axis(side = 1, at = c(min(x.ticks), max(info$date)), lwd = 1, lwd.ticks = 0, labels = F)
# axis(side = 1, at = x.ticks[seq.int(1,20,2)], labels = F, tck = -.025, lwd = 0, lwd.ticks = 1)
# axis(side = 1, at = x.ticks[seq.int(1,20,2)], labels = year(x.ticks)[seq.int(1,20,2)], lwd = 0, las = 1, line = -.5)
# axis(side = 2, labels = F, line = -.75, tck = -.025)
# axis(side = 2, lwd = F, las = 2, line = -1.25)
# mtext(text = "Relative Abundance (%)", side = 2, line = 2)
# dev.off()
# # ----

# acI-B nucleotide diversity ----

par(mar = c(1.2,3.27,.3,0))
plot(x = info$date, y = info$nucl_diversity, type = "n", ann = F, axes = F)
fill.under.lines(X = info$date, Y = info$nucl_diversity, YAxisMin = 0, Color = adjustcolor(col = "#a9bfdb", alpha.f = 1))
lines(x = info$date, y = info$nucl_diversity, col = "#314250")
points(x = info$date, y = info$nucl_diversity, pch = 16, col = adjustcolor("#314250", alpha.f = .7), cex = .5)
axis(side = 1, at = c(min(x.ticks), max(info$date)), lwd = 1, lwd.ticks = 0, labels = F)
axis(side = 1, at = x.ticks, labels = F, tck = -.05, lwd = 0, lwd.ticks = 1)
axis(side = 1, at = x.ticks[1], labels = year(x.ticks)[1], lwd = 0, las = 1, line = -.75, hadj = .25)
axis(side = 1, at = x.ticks[2], labels = year(x.ticks)[2], lwd = 0, las = 1, line = -.75, hadj = .5)
axis(side = 1, at = x.ticks[3], labels = year(x.ticks)[3], lwd = 0, las = 1, line = -.75, hadj = .75)
axis(side = 2, at = c(.008,.013,.018), labels = F, line = -.25, tck = -.05)
axis(side = 2, at = c(.008,.013,.018), labels = c(".008",".013",".018"), lwd = 0, las = 2, line = -.75)
mtext(text = expression(paste("Nucl. Div. (",pi,")")), side = 2, line = 2.1)

plot(x = info$date, y = info$nucl_diversity, type = "n", ann = F, axes = F, xlim = parse_date_time(x = c(2011,2014), orders = "y"), ylim = c(.011,.017))
fill.under.lines(X = info$date, Y = info$nucl_diversity, YAxisMin = 0, Color = adjustcolor(col = "#a9bfdb", alpha.f = 1))
lines(x = info$date, y = info$nucl_diversity, col = "#314250")
points(x = info$date, y = info$nucl_diversity, pch = 16, col = adjustcolor("#314250", alpha.f = .7), cex = .5)
text(x = info$date, y = info$nucl_diversity, labels = info$date, cex = .6, srt = 30)

# ----


# RRR
# Make NMDS plots with base R, force them to be square etc
# make sure lines are actually going in order!
# Use Euclidean distance because total SNVs does not have a max carrying capacity, and increases after 2012

library(data.table)
library(lubridate)
library(viridisLite)

per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-21_3300043464_group3_bin69_SNVs.tsv.gz" # B
genome.info <- fread(input = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
genome.info[genome ==  "ME2011-09-21_3300043464_group3_bin69",above.breadth.cutoff := (breadth / breadth_expected) >= .5, ] # B
num.threads <- 1

# ---- 

x <- "character"
names(x) <- "date"
snv.instrain <- fread(file = per.genome.snv.file, nThread = num.threads, colClasses = x)
my.genome <- snv.instrain[1, genome]

snv.instrain

genome.info <- genome.info[above.breadth.cutoff == TRUE & coverage_median > 10]

snv.instrain <- snv.instrain[sample %in% genome.info$sample]

snvs <- dcast(data = snv.instrain, formula = scaffold + position + gene + mutation_type ~ sample, value.var = "ref_freq", fill = 1)

snvs[ ,snv.num := paste0("snv",1:nrow(snvs))]
snv.mat <- as.matrix(snvs[ ,-c(1:4)], rownames = "snv.num")
snv.mat <- t(snv.mat)

snv.dist.euc <- vegdist(x = snv.mat, method = "euclidean")
snv.nmds.euc <- metaMDS(comm = snv.dist.euc, trymax = 100)
stressplot(snv.nmds.euc)

dates.key <- unique(snv.instrain[ ,.(sample, date, year, yday, season, invasion)])

my.nmds <- data.table("sample" = rownames(snv.nmds.euc$points),"x" = snv.nmds.euc$points[ ,1], "y" = snv.nmds.euc$points[ ,2])

my.nmds <- merge(x = my.nmds, y = dates.key, by = "sample")
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
text(x = my.nmds$x[index.2012], y = my.nmds$y[index.2012], labels = my.nmds$date[index.2012], cex = .7)

par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 2, to = 2, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)

