# RRR
# the NMDS objects have already been generated for each genome by 1-calculate_distance_matrix_and_NMDS.R
# that script subset to >10 median coverate already, so low abund samples are already missing from it
# it also made versions with all SNVs and with nonsynonymous-only SNVs
# take the NMDS object, subset the sample key to match it
# make NMDS plots colored by:
#   - pre/post 2012
#   -long-term gradual change
#   -amount of positive selection
#   -amount of negative selection
#   -genome abundance
#   -genome nucleotide diversity
#   -season

# NOTE: make a bunch of output folders before running!

# ---- set-up ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(viridis))

userprefs <- commandArgs(trailingOnly = TRUE)
nmds.file <- userprefs[1]
sample.key.file <- userprefs[2]
selection.summary.file <- userprefs[3]
genome.file <- userprefs[4]
tax.file <- userprefs[5]
output.folder.tables <- userprefs[6]
output.folder.plots <- userprefs[7]
threads <- as.numeric(userprefs[8])

# # example paths for local testing
# nmds.file <- "data/2023-11-01_multidimensional_SNV_analysis/example_input_files/ME2015-06-10_3300042356_group6_bin86_all_SNV_euclidean_distance_NMDS_object.rds"
# sample.key.file <- "data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv"
# selection.summary.file <- "data/2023-10-26_selection_per_sample_summaries/example_data/ME2015-06-10_3300042356_group6_bin86_selection_summary.tsv.gz"
# genome.file <- "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz"
# tax.file <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"
# output.folder.tables <- "data/2023-11-01_multidimensional_SNV_analysis/NMDS_example_plots"
# output.folder.plots <- "data/2023-11-01_multidimensional_SNV_analysis/NMDS_example_plots"
# threads <- 1


# ---- import and add color scales ----

my.genome <- sub(pattern = "^.*ME", replacement = "ME", x = selection.summary.file)
my.genome <- sub(pattern = "_selection_summary.*$", replacement = "", x = my.genome)
cat("Processing", my.genome,"\n")

# not all genomes HAVE an NMDS object, they had to have enough coverage. skip if missing the file:
if (!file.exists(nmds.file)){
  quit(save = "no", status = 0)
}
# not all genomes HAVE a selection summary file, they had to have enough coverage and I think I removed the GD samples. skip if missing the file:
if (!file.exists(selection.summary.file)){
  quit(save = "no", status = 0)
}

nmds.obj <- readRDS(nmds.file)

sample.key <- fread(file = sample.key.file, nThread = 1, colClasses = c("character","character","numeric","numeric","character","character"))
sample.key[ ,date := parse_date_time(x = date, orders = "ymd")]

# date colors spaced by true time, not number of samples
sample.key$julian.date <- as.numeric(julian(sample.key$date))
color.key.date <- data.table("julian.date" = min(sample.key$julian.date):max(sample.key$julian.date))
color.key.date[ ,color.date := viridis(n = nrow(color.key.date), option = "viridis")]
sample.key <- merge(x = sample.key, y = color.key.date, by = "julian.date", all.x = TRUE, all.y = FALSE)
color.key.date[ ,date := as.Date(x = julian.date, origin = as.Date("1970-01-01"))]
color.key.date[ ,year := year(date)]

color.key.season <- data.table("season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), "color.season" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))
sample.key <- merge(x = sample.key, y = color.key.season, by = "season")
color.key.season[ ,season := sub(pattern = "Ice.On", replacement = "Ice-On", x = season)]
color.key.season[ ,season := sub(pattern = "\\.", replacement = " ", x = season)]

color.key.2012 <- data.table("year" = 2000:2019, "color.2012" = c(rep("#b89996", 11), "purple","red","tan2","dodgerblue4", rep("#7592b2", 5)))
sample.key <- merge(x = sample.key, y = color.key.2012, by = "year")

selection.summary <- fread(file = selection.summary.file, nThread = threads, colClasses = c("character","character","character","character","numeric","numeric","numeric"))
selection.summary <- selection.summary[ ,c("sample","Npos","Nneg")]

if (nrow(selection.summary) > 0){
  color.key.pos <- data.table("Npos" = min(selection.summary$Npos):max(selection.summary$Npos))
  color.key.pos[ ,color.pos := viridis(n = nrow(color.key.pos), option = "rocket")]
  selection.summary <- merge(x = selection.summary, y = color.key.pos, by = "Npos", all.x = TRUE, all.y = FALSE)
  
  color.key.neg <- data.table("Nneg" = min(selection.summary$Nneg):max(selection.summary$Nneg))
  color.key.neg[ ,color.neg := viridis(n = nrow(color.key.neg), option = "cividis")]
  selection.summary <- merge(x = selection.summary, y = color.key.neg, by = "Nneg", all.x = TRUE, all.y = FALSE)
  
  sample.key <- merge(x = sample.key, y = selection.summary, by = "sample", all = TRUE)
}

tax <- readRDS(file = tax.file)
tax <- tax[ ,c("bin.full.name","completeness","contamination","num.in.cluster","domain","phylum","class","order","family","genus","species")]
colnames(tax)[1] <- "genome"
tax <- as.data.table(tax)
tax <- tax[genome == my.genome]
tax.label <- paste(tax[ ,.(phylum,class,order,family,genus,species)], collapse = " ")

nmds.table <- data.table("genome" = my.genome, "sample" = row.names(nmds.obj$points),"x" = nmds.obj$points[ ,1], "y" = nmds.obj$points[ ,2])

nmds.table <- merge(x = nmds.table, y = sample.key, by = c("sample"), all.x = TRUE, all.y = FALSE)

my.lims <- c(-max(abs(c(nmds.table$x, nmds.table$y)))-5, max(abs(c(nmds.table$x, nmds.table$y)))+5)

# ---- make and export stressplot ----

pdf(file = file.path(output.folder.plots, paste0(tax.label," ",my.genome," 1-stressplot.pdf")), width = 6.5, height = 6.5)
stressplot(nmds.obj)
mtext(text = tax.label, side = 3, adj = 0, cex = .5, line = 3)
mtext(text = my.genome, side = 3, adj = 0, cex = 1, line = 1.5)
dev.off()

# ---- plot gradual change ----

pdf(file = file.path(output.folder.plots, paste0(tax.label," ",my.genome," 2-gradual_change.pdf")), width = 4.2, height = 3)

par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = nmds.table$x, y = nmds.table$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = nmds.table$x, y = nmds.table$y, pch = 16, col = adjustcolor(nmds.table$color.date, alpha.f = .7), cex = .7)

par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(2.5,.1,2.5,5))
plot(x = rep(1,length(color.key.date$julian.date)), y = color.key.date$julian.date, xlim = c(0,1), col = color.key.date$color.date, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = color.key.date$julian.date, y1 = color.key.date$julian.date, col = color.key.date$color.date, lwd = 5)
index.labels <- which(yday(color.key.date$date) == 1)
index.labels <- index.labels[seq.int(1,length(index.labels),2)]
axis(side = 4, at = color.key.date$julian.date[index.labels], labels = color.key.date$year[index.labels], xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = color.key.date$julian.date[index.labels], xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
mtext(text = "Sample Date", adj = 0, line = 1, cex = 1, at = 0)

dev.off()


pdf(file = file.path(output.folder.plots, paste0(tax.label," ",my.genome," 3-gradual_change-lines.pdf")), width = 4.2, height = 3)

par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = nmds.table$x, y = nmds.table$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = nmds.table$x, y = nmds.table$y, pch = 16, col = adjustcolor(nmds.table$color.date, alpha.f = .7), cex = .7)
lines(x = nmds.table$x, y = nmds.table$y, pch = 16, col = adjustcolor("grey40", alpha.f = .5), lwd = .5)

par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(2.5,.1,2.5,5))
plot(x = rep(1,length(color.key.date$julian.date)), y = color.key.date$julian.date, xlim = c(0,1), col = color.key.date$color.date, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = color.key.date$julian.date, y1 = color.key.date$julian.date, col = color.key.date$color.date, lwd = 5)
index.labels <- which(yday(color.key.date$date) == 1)
index.labels <- index.labels[seq.int(1,length(index.labels),2)]
axis(side = 4, at = color.key.date$julian.date[index.labels], labels = color.key.date$year[index.labels], xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = color.key.date$julian.date[index.labels], xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
mtext(text = "Sample Date", adj = 0, line = 1, cex = 1, at = 0)

dev.off()

# ---- plot 2012 ----

pdf(file = file.path(output.folder.plots, paste0(tax.label," ",my.genome," 4-2012.pdf")), width = 4.2, height = 3)

par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = nmds.table$x, y = nmds.table$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = nmds.table$x, y = nmds.table$y, pch = 16, col = adjustcolor(nmds.table$color.2012, alpha.f = .7), cex = .7)

par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(2.5,.1,2.5,5))
par(mar = c(.1,.1,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 0, to = 0, along.with = 1:6), 
       col = c(adjustcolor(col = c("#b89996", "purple"), alpha.f = .7), "red", adjustcolor(col = c("tan2","dodgerblue4","#7592b2"), alpha.f = .7)), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:6), x = seq(from = 1.75, to = 1.75, along.with = 1:6),
     labels = c("2000 - 2010", "2011","2012","2013","2014","2015 - 2019"), adj = 0, xpd = NA)

dev.off()

# ---- plot seasons ----

pdf(file = file.path(output.folder.plots, paste0(tax.label," ",my.genome," 5-seasons.pdf")), width = 4.2, height = 3)

par(fig = c(0,3/4.2,0,1))
par(mar = c(1.5,1.5,1.5,1.5))
plot(x = nmds.table$x, y = nmds.table$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
box()
mtext(text = "NMDS Axis 1", side = 1, line = .25)
mtext(text = "NMDS Axis 2", side = 2, line = .25)
points(x = nmds.table$x, y = nmds.table$y, pch = 16, col = adjustcolor(nmds.table$color.season, alpha.f = .7), cex = .7)

par(fig = c(3/4.2,1,0,1), new = T)
par(mar = c(2.5,.1,2.5,5))
par(mar = c(.1,.075,.1,.1))
plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)
points(y = seq(from = 8, to = 4, along.with = 1:nrow(color.key.season)), x = seq(from = 0, to = 0, along.with = 1:nrow(color.key.season)), 
       col = adjustcolor(col = color.key.season$color.season, alpha.f = .7), 
       pch = 16, cex = 2, xpd = NA)
text(y = seq(from = 8, to = 4, along.with = 1:nrow(color.key.season)), x = seq(from = 1.25, to = 1.25, along.with = 1:nrow(color.key.season)),
     labels = color.key.season$season, adj = 0, xpd = NA)

dev.off()

# ---- plot selection ----

if (nrow(selection.summary) > 0){
  
  pdf(file = file.path(output.folder.plots, paste0(tax.label," ",my.genome," 6-positive_selection.pdf")), width = 4.2, height = 3)
  
  par(fig = c(0,3/4.2,0,1))
  par(mar = c(1.5,1.5,1.5,1.5))
  plot(x = nmds.table$x, y = nmds.table$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
  box()
  mtext(text = "NMDS Axis 1", side = 1, line = .25)
  mtext(text = "NMDS Axis 2", side = 2, line = .25)
  points(x = nmds.table$x, y = nmds.table$y, pch = 16, col = adjustcolor(nmds.table$color.pos, alpha.f = .7), cex = .7)
  
  par(fig = c(3/4.2,1,0,1), new = T)
  par(mar = c(5,0,5,5))
  plot(x = rep(1,length(color.key.pos$Npos)), y = color.key.pos$Npos, xlim = c(0,1), col = color.key.pos$color.pos, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
  segments(x0 = 0, x1 = 1, y0 = color.key.pos$Npos, y1 = color.key.pos$Npos, col = color.key.pos$color.pos, lwd = 10)
  label.locs <- pretty(x = color.key.pos$Npos)
  label.locs <- label.locs[label.locs <= max(color.key.pos$Npos) & label.locs >= min(color.key.pos$Npos)]
  label.locs <- round(label.locs, digits = 0)
  label.locs <- unique(label.locs)
  axis(side = 4, at = label.locs, xpd = T, lwd = 0, las = 2, line = -.5, cex.axis = 1)
  axis(side = 4, at = label.locs, xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
  mtext(text = "Genes under\npositive selection", adj = 0, line = 1, cex = 1, at = -.5)
  
  dev.off()
  
  
  pdf(file = file.path(output.folder.plots, paste0(tax.label," ",my.genome," 7-negative_selection.pdf")), width = 4.2, height = 3)
  
  par(fig = c(0,3/4.2,0,1))
  par(mar = c(1.5,1.5,1.5,1.5))
  plot(x = nmds.table$x, y = nmds.table$y, xlim = my.lims, ylim = my.lims, xaxs = "i", yaxs = "i", ann = F, axes = F, type = "n")
  box()
  mtext(text = "NMDS Axis 1", side = 1, line = .25)
  mtext(text = "NMDS Axis 2", side = 2, line = .25)
  points(x = nmds.table$x, y = nmds.table$y, pch = 16, col = adjustcolor(nmds.table$color.neg, alpha.f = .7), cex = .7)
  
  par(fig = c(3/4.2,1,0,1), new = T)
  par(mar = c(5,0,5,5))
  plot(x = rep(1,length(color.key.neg$Nneg)), y = color.key.neg$Nneg, xlim = c(0,1), col = color.key.neg$color.neg, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
  segments(x0 = 0, x1 = 1, y0 = color.key.neg$Nneg, y1 = color.key.neg$Nneg, col = color.key.neg$color.neg, lwd = 10)
  label.locs <- pretty(x = color.key.neg$Nneg)
  label.locs <- label.locs[label.locs <= max(color.key.neg$Nneg) & label.locs >= min(color.key.neg$Nneg)]
  label.locs <- round(label.locs, digits = 0)
  label.locs <- unique(label.locs)
  axis(side = 4, at = label.locs, xpd = T, lwd = 0, las = 2, line = -.5, cex.axis = 1)
  axis(side = 4, at = label.locs, xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
  mtext(text = "Genes under\npurifying selection", adj = 0, line = 1, cex = 1, at = -.5)
  
  dev.off()
  
}

# ---- export plot table ----

fwrite(x = nmds.table[ ,date := as.character(date)], file = file.path(output.folder.tables, paste0(my.genome,"_NMDS_table.tsv.gz")), nThread = threads, sep = "\t")

# ---- end ----
