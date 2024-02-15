# RRR
# Make the abundance over time plot a single color

library(data.table)
library(lubridate)

coverm <- fread(input = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
genome.info <- fread(input = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")

my.genome <- "ME2011-09-21_3300043464_group3_bin69" # acI-B
# my.genome <- "ME2016-07-20_3300033996_group7_bin32" # acI-C
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
pdf(file = "figures/2023-08-10_UCSB_plots/nuc_div_acI-B_long-term.pdf",width = 3.1, height = 1.36)
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
dev.off()
# ----
