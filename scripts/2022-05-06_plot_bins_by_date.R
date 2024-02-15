# RRR
# can't log on to IMG, but will make the script work with an old version of the data


library(readxl)
library(lubridate)
source("pop/scripts-generic_plots/plot_sample_dates_fancy.R")
source("pop/scripts-processing/5C_subset_by_taxon.R")

limony <- readRDS("pop/data/limony-IRD/2021-08-25_processing/5A_taxlist_samples_6000.rds")
bigtyme <- read_excel(path = "data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx")
bin.stats <- read_excel(path = "data/2022-01-05_IMG_IDs_for_Actino_bins/IMG_bins_Actinobacteria.xlsx")
colnames(bin.stats)
colnames(bigtyme)

my.bins <- merge(x = bin.stats, y = bigtyme, by.x = "IMG Genome ID", by.y = "taxon_oid", all.x = TRUE, all.y = FALSE)
my.bins$Date <- parse_date_time(x = paste(my.bins$Year,my.bins$Month, my.bins$Day), orders = "ymd", tz = "Etc/GMT-5")

# subset to taxon of interest -bins
index <- grep(pattern = "Nanopelagicales", x = my.bins$`GTDBTK Lineage`, ignore.case = T, value = F)
my.data <- my.bins[index, ]

# subset to taxon of interest -limony
lim <- flatten.to.single.level(my.list = limony, finest.level = "acI")
index <- get.taxon.indeces(my.list = lim, taxa = "acI")
lim <- subset.by.taxa(my.list = lim, keep.index = index, renormalize = F)
my.tags <- lim$av[1, ]
my.dates <- convert.colnames.to.dates(sample.names = names(my.tags))

# look at number per day:
my.data$Count <- 1
my.data <- aggregate(x = my.data$Count, by = list(my.data$Date), FUN = sum)
colnames(my.data) <- c("Date","Num.Bins")

my.cex <- get.scaled.point.cex(min.cex = .5, max.cex = 3, my.values = my.data$Num.Bins)

plot.info <- set.up.plot(dates.vector = my.bins$Date, year.lines = T, month.lines = T)
add.dates(date.vector = my.bins$Date, year.labs = plot.info, point.type = 21, point.col = "white", line.col = "black", point.cex = .5)
add.dates(date.vector = my.data$Date, year.labs = plot.info, point.type = 21, point.col = "red", line.col = "black", point.cex = my.cex)
mtext(text = "All Nanopelagicales bin in red, size ~ number, max size = 13 bins")

# overlay acI abundance

tags.offset <- get.year.scaling.and.offset(my.max = max(my.tags), my.min = min(my.tags))
my.tags <- get.new.x.y(x = my.dates, y = my.tags, axis.year = 2000, year.offset = tags.offset$offset, year.scaling = tags.offset$scaling)
for (y in unique(my.tags$year)){
  i <- which(my.tags$year == y)
  lines(x = my.tags$x[i], y = my.tags$y[i])
}

# filter bins more by completeness

index <- grep(pattern = "Nanopelagicales", x = my.bins$`GTDBTK Lineage`, ignore.case = T, value = F)
my.data <- my.bins[index, ]
hist(my.data$`Bin Completeness`)

index <- which(my.data$`Bin Completeness` >= 70)
my.data <- my.data[index, ]

# look at number per day:
my.data$Count <- 1
my.data <- aggregate(x = my.data$Count, by = list(my.data$Date), FUN = sum)
colnames(my.data) <- c("Date","Num.Bins")

max(my.data$Num.Bins)
my.cex <- get.scaled.point.cex(min.cex = 1, max.cex = 2, my.values = my.data$Num.Bins)

plot.info <- set.up.plot(dates.vector = my.bins$Date, year.lines = T, month.lines = T)
add.dates(date.vector = my.bins$Date, year.labs = plot.info, point.type = 21, point.col = "white", line.col = "black", point.cex = .5)
add.dates(date.vector = my.data$Date, year.labs = plot.info, point.type = 21, point.col = "red", line.col = "black", point.cex = my.cex)
mtext(text = "All Nanopelagicales  >= 90% completeness in red, size ~ number, max size = 1 bins")

# overlay acI abundance
for (y in unique(my.tags$year)){
  i <- which(my.tags$year == y)
  lines(x = my.tags$x[i], y = my.tags$y[i])
}

# gather some quick number of bin stats
index <- grep(pattern = "Nanopelagicales", x = my.bins$`GTDBTK Lineage`, ignore.case = T, value = F)
my.data <- my.bins[index, ]
hist(my.data$`Bin Completeness`)

length(which(my.data$`Bin Completeness` >= 50)) # 599 that were MQ and completeness >= 50
length(which(my.data$`Bin Completeness` >= 60)) # 291
length(which(my.data$`Bin Completeness` >= 70)) # 109
length(which(my.data$`Bin Completeness` >= 80)) # 25
length(which(my.data$`Bin Completeness` >= 90)) # 3
length(which(my.data$`Bin Completeness` >= 95)) # 1



