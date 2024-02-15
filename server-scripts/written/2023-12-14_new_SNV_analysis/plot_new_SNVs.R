# RRR
# make some simple plots of total SNVs and new SNV for all the abundant genomes.
# just to look through, probably will not be included in the paper.
# worry about classifying immigration events later.

# ---- set-up ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(lubridate))

userprefs <- commandArgs(trailingOnly = TRUE)
per.genome.stats.file <- userprefs[1]
tax.file <- userprefs[2]
num.threads <- as.numeric(userprefs[3])
output.folder <- userprefs[4]

# # local paths testing
# per.genome.stats.file <- "data/2023-12-14_new_SNV_testing/ME2014-08-28_3300035662_group5_bin212_new_and_total_SNVs.tsv.gz"
# tax.file <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"
# num.threads <- 1
# output.folder <- "data/2023-12-14_new_SNV_testing/"

# ---- prepare data ----

if (!file.exists(per.genome.stats.file)){
  quit(save = "no", status = 0)
}

snv.stats <- fread(file = per.genome.stats.file, colClasses = c("Date" = "character"), nThread = num.threads)
snv.stats[ ,Date := parse_date_time(x = Date, orders = "ymd")]

my.genome <- snv.stats$genome[1]
cat("Processing",my.genome,"\n")

tax <- readRDS(file = tax.file)
tax <- as.data.table(tax)
tax <- tax[bin.full.name == my.genome]
tax.label <- paste(tax$phylum,tax$class,tax$order,tax$family,tax$genus,tax$species)

# ---- for plots ----

x.ax.labs <- seq.int(from = min(snv.stats$year), to = max(snv.stats$year), by = 2)
x.ax.locs <- parse_date_time(x = paste(x.ax.labs,"-1-1"), orders = "ymd")

# exclude the high initial values
q3 <- quantile(x = snv.stats$New, na.rm = T)  # get the high values
high.ones <- snv.stats[New > q3[4], ]
start.date <- high.ones$Date[cumsum(high.ones$Date != snv.stats$Date[-1][1:length(high.ones$Date)]) == 1]
if (length(start.date) == 0){
  start.date <- snv.stats$Date[-1][length(high.ones$Date) + 1]
}
zoom.index <- snv.stats$Date >= start.date

# and find any high new-SNV outliers after it levels off
my.outliers <- boxplot(x = snv.stats[zoom.index ,New], plot = F)
my.outliers <- my.outliers$out[my.outliers$out > my.outliers$stats[3,1]] # onlly the high outliers


# ---- make plot ----

pdf(file = file.path(output.folder,paste(tax.label,my.genome,".pdf")), width = 8, height = 5)

par(fig = c(0,1,.7,1), mar = c(.1,5,3,.5))
plot(x = snv.stats$Date, y = snv.stats$tot.SNVs, type = "l", ann = F, axes = F)
box()
axis(side = 2, labels = F, lwd = 0, lwd.ticks = 1)
axis(side = 2, labels = T, lwd = 0, las = 2, line = -.25)
mtext(side = 2, text = "Tot. SNVs", line = 4, adj = .5)
mtext(text = tax.label, side = 3, line = 1.75, adj = 0, cex = .7)
mtext(text = my.genome, side = 3, line = .75, adj = 0)

par(fig = c(0,1,.4,.7), mar = c(3,5,.1,.5), new = T)
plot(x = snv.stats$Date, y = snv.stats$New, type = "l", ann = F, axes = F)
box()
axis(side = 1, at = x.ax.locs, labels = F, lwd = 0, lwd.ticks = 1)
axis(side = 1, at = x.ax.locs, labels = x.ax.labs, lwd = 0, line = -.75)
axis(side = 2, labels = F, lwd = 0, lwd.ticks = 1)
axis(side = 2, labels = T, lwd = 0, las = 2, line = -.25)
mtext(side = 2, text = "New SNVs", line = 4, adj = -.5)

par(fig = c(0,1,0,.4), mar = c(2.5,5,.1,.5), new = T)
plot(x = snv.stats[zoom.index, Date], y = snv.stats[zoom.index ,New], type = "l", ann = F, axes = F)
points(x = snv.stats[New %in% my.outliers, Date], y = snv.stats[New %in% my.outliers, New], pch = 19, col = "red")
box(col = "red")
axis(side = 1, at = x.ax.locs, labels = F, lwd = 0, lwd.ticks = 1, col = "red")
axis(side = 1, at = x.ax.locs, labels = x.ax.labs, lwd = 0, line = -.75, col.axis = "red")
axis(side = 2, labels = F, lwd = 0, lwd.ticks = 1)
axis(side = 2, labels = T, lwd = 0, las = 2, line = -.25)
mtext(side = 2, text = "New SNVs", line = 4)
mtext(side = 2, text = "(zoomed)", line = 3, col = "red")
points(x = parse_date_time(paste(year(min(snv.stats[zoom.index, Date])),"-1-1"),orders = "ymd"), y = max(snv.stats[zoom.index ,New]) + max(snv.stats[zoom.index ,New]) / 7, col = "red", xpd = NA, pch = 19)
text(x = parse_date_time(paste(year(min(snv.stats[zoom.index, Date])),"-6-1"),orders = "ymd"), y = max(snv.stats[zoom.index ,New]) + max(snv.stats[zoom.index ,New]) / 7, labels = "outliers (immigration events?)", xpd = NA, adj = 0, col = "red")
dev.off()

# ----


