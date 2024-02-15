# RRR
# Make the abundance over time plot a single color

library(data.table)
library(lubridate)

coverm <- fread(input = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
genome.info <- fread(input = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")

my.genome <- "ME2011-09-21_3300043464_group3_bin69" # acI-B
my.genome <- "ME2016-07-20_3300033996_group7_bin32" # acI-C
my.genome <- "ME2011-09-04_3300044729_group3_bin142" # acI-A
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

x.ticks <- parse_date_time(x = 2000:2019, orders = "y")

# acI-B abund ----
pdf(file = "figures/2023-08-05_ESA_plots/abund_acI-B_long-term.pdf",width = 6.5, height = 3)
par(mar = c(1.75,3,1,0))
plot(x = info$date, y = info$rel.abund, type = "n", ann = F, axes = F)
fill.under.lines(X = info$date, Y = info$rel.abund, YAxisMin = 0, Color = adjustcolor(col = "#b89996", alpha.f = .5))
lines(x = info$date, y = info$rel.abund, col = "#b89996")
points(x = info$date, y = info$rel.abund, pch = 16, col = adjustcolor("#b89996", alpha.f = .7))
axis(side = 1, at = c(min(x.ticks), max(info$date)), lwd = 1, lwd.ticks = 0, labels = F)
axis(side = 1, at = x.ticks[seq.int(1,20,2)], labels = F, tck = -.025, lwd = 0, lwd.ticks = 1)
axis(side = 1, at = x.ticks[seq.int(1,20,2)], labels = year(x.ticks)[seq.int(1,20,2)], lwd = 0, las = 1, line = -.5)
axis(side = 2, labels = F, line = -.75, tck = -.025)
axis(side = 2, lwd = F, las = 2, line = -1.25)
mtext(text = "Relative Abundance (%)", side = 2, line = 2)
dev.off()
# ----

# acI-B nucleotide diversity ----
pdf(file = "figures/2023-08-05_ESA_plots/nuc_div_acI-B_long-term.pdf",width = 6.5, height = 3)
par(mar = c(1.75,3.5,1,0))
plot(x = info$date, y = info$nucl_diversity, type = "n", ann = F, axes = F)
fill.under.lines(X = info$date, Y = info$nucl_diversity, YAxisMin = 0, Color = adjustcolor(col = "#7592b2", alpha.f = .5))
lines(x = info$date, y = info$nucl_diversity, col = "#7592b2")
points(x = info$date, y = info$nucl_diversity, pch = 16, col = adjustcolor("#7592b2", alpha.f = .7))
axis(side = 1, at = c(min(x.ticks), max(info$date)), lwd = 1, lwd.ticks = 0, labels = F)
axis(side = 1, at = x.ticks[seq.int(1,20,2)], labels = F, tck = -.025, lwd = 0, lwd.ticks = 1)
axis(side = 1, at = x.ticks[seq.int(1,20,2)], labels = year(x.ticks)[seq.int(1,20,2)], lwd = 0, las = 1, line = -.5)
axis(side = 2, labels = F, line = -.75, tck = -.025)
axis(side = 2, lwd = F, las = 2, line = -1.25)
mtext(text = expression(paste("Nucleotide Diversity (",pi,")")), side = 2, line = 2.3)
dev.off()
# ----

# ---- plot set-up ggplot ----
library(ggplot2)
library(patchwork)

info[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]

p.abund <- ggplot(data = info, aes(x = yday, y = rel.abund))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .2)+
  geom_smooth(color = "black", se = T)+
  geom_point(aes(color = season, group = year), size = 1.5, alpha = .8)+
  geom_smooth(color = "black", se = F)+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))+
  scale_y_continuous(name = "Relative\nAbundance(%)")+
  guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5)))


p.nuc <- ggplot(data = info, aes(x = yday, y = nucl_diversity))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .2)+
  geom_smooth(color = "black", se = T)+
  geom_point(aes(color = season, group = year), size = 1.5, alpha = .8)+
  geom_smooth(color = "black", se = F)+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))+
  scale_y_continuous(name = expression(paste("Nucleotide\nDiversity (",pi,")")))+
  guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5)))

pdf(file = "figures/2023-08-05_ESA_plots/abund_acI-C_seasons.pdf",width = 6.5, height = 3)
pdf(file = "figures/2023-08-05_ESA_plots/abund_acI-A_seasons.pdf",width = 6.5, height = 3)
p.abund
dev.off()

pdf(file = "figures/2023-08-05_ESA_plots/nuc_div_abund_acI-C_seasons.pdf",width = 6.5, height = 3)
pdf(file = "figures/2023-08-05_ESA_plots/nuc_div_abund_acI-A_seasons.pdf",width = 6.5, height = 3)
p.abund / p.nuc + plot_layout(guides = "collect")
dev.off()
