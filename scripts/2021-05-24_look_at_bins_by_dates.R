# Look at bin rank abunds.... maybe I should be doing this in pop?

library(lubridate)

bins <- readRDS(file = "data/2021-05-22_bin_stats/all_img_bins.rds")
asms <- readRDS(file = "data/2021-05-24_metaG_stats/all_img_assemblies.rds")

colnames(bins)
unique(bins$IMG.Phylum)
unique(bins$GTDB.Phylum)

colnames(asms)

# ---- functions ----

plot.bin.rank.abund <- function(bins, tax.col = 17, nbar = NA){
  hq <- bins[bins$Bin.Quality == "HQ", ]
  mq <- bins[bins$Bin.Quality == "MQ", ]
  
  count.vect <- rep.int(1,nrow(hq))
  hq <- aggregate(x = count.vect, by = list(hq[ ,tax.col]), FUN = "sum")
  names(hq) <- c("tax","HQ.bins")
  
  count.vect <- rep.int(1,nrow(mq))
  mq <- aggregate(x = count.vect, by = list(mq[ ,tax.col]), FUN = "sum")
  names(mq) <- c("tax","MQ.bins")
  
  bin.mat <- merge(x = hq, y = mq, by = "tax", all = TRUE)
  row.names(bin.mat) <- bin.mat$tax
  bin.mat <- bin.mat[ ,2:3]
  bin.mat <- as.matrix(bin.mat)
  index <- is.na(bin.mat)
  bin.mat[index] <- 0
  index <- order(rowSums(bin.mat), decreasing = T)
  bin.mat <- bin.mat[index, ]
  bin.mat <- t(bin.mat)
  tot.names <- ncol(bin.mat)
  if (!is.na(nbar) & nbar < tot.names){
    bin.mat <- bin.mat[ ,1:nbar]
  }
  barplot(height = bin.mat, names.arg = colnames(bin.mat), beside = F, las = 2, 
          col = c("thistle3","thistle1"), legend = c("HQ","MQ"), cex.names = .9)
  mtext(text = paste0("Rank Abundance of Total Bins by ", colnames(bins)[tax.col],
                      "\n(",nrow(bins), " total bins, ",tot.names," total ",colnames(bins)[tax.col] ,")"), line = -1.5)
}


# From pop:
convert.colnames.to.dates <- function(sample.names){
  sample.names <- substr(x = sample.names, start = 3, stop = 11)
  sample.names <- parse_date_time(x = sample.names, orders = "dmy", tz = "Etc/GMT-5")
  return(sample.names)
}

get.plot.x.vals <- function(my.y, label.with, step.size = 2, label.btwn.ticks = T, force.x.min.YY.MM.DD = FALSE, force.x.max.YY.MM.DD = FALSE){
  # label.with can be "year", "month", or "day"
  # label location is in middle of range, tick loc is at extremes. ie label is at June 15, ticks are at June 1 and June 31
  # my.y is the flat.list structure, subset to n rows, that also includes ymax
  
  # ---- functions ----
  
  check.function.input.format <- function(label.with, step.size){
    if (label.with != "year" & label.with != "month" & label.with != "day"){
      return("label.with should be exactly either \"year\", \"month\", or \"day\"")
    }
    if (all.equal(floor(step.size), step.size) != TRUE){
      return(cat("step size has to be an integer. work around examples:",
                 "for half-year ticks, add another axis of ticks that is transformed by date_decimal(x.ticks + .5)",
                 "for ticks inside of months, repeat this function call to get a vector of day-ticks with desired step-size", sep = "\n"))
    }
    return(cat("Creating x-axis labeled at every", step.size, label.with, "\n"))
  }
  
  get.tick.locs <- function(num.ticks, label.with, step.size, pts.span){
    x.ticks <- seq.int(from = 0, length.out = num.ticks, by = step.size)
    x.ticks <- lapply(X = as.list(x.ticks), FUN = period, units = label.with) # why lapply and the do.call?  https://stackoverflow.com/questions/15659783/why-does-unlist-kill-dates-in-r
    x.ticks <- do.call(what = "c", args = x.ticks)                             # b/c need apply b/c period sums vector arguments, need lapply+do.call b/c sapply calls unlist which drops date attributes
    x.ticks <- int_start(pts.span) %m+% x.ticks
    return(x.ticks)
  }
  
  get.x.labs <- function(label.with, x.ticks, offset.labs){
    if (label.with == "month"){
      x.labs <- month(x.ticks, label = TRUE, abbr = TRUE) # specify want text not numbers
    }else{
      x.labs <- sapply(X = x.ticks, FUN = label.with)
    }
    
    if (offset.labs){
      x.labs <- x.labs[-(length(x.labs))]
    }
    
    return(x.labs)
  }
  
  get.lab.locs <- function(label.with, x.ticks, offset.labs){
    if (offset.labs){
      if (label.with == "year"){
        x.lab.locs <- date_decimal(x.labs + .5)
      }else if(label.with == "month"){
        x.lab.locs <- x.ticks[-num.ticks] %m+% days(14)
      }else if(label.with == "day"){
        x.lab.locs <- x.ticks[-num.ticks] %m+% hours(12)
      }
    }else{
      x.lab.locs <- x.ticks
    }
    return(x.lab.locs)
  }
  
  
  # ---- actions ----
  
  check.function.input.format(label.with = label.with, step.size = step.size)
  
  dates.vector <- convert.colnames.to.dates(sample.names = colnames(my.y$bq))
  
  if (force.x.min.YY.MM.DD == FALSE){
    x.min <- min(dates.vector)
  }else{
    x.min <- parse_date_time(x = force.x.min.YY.MM.DD, orders = "ymd", tz = "Etc/GMT-5")
  }
  if (force.x.max.YY.MM.DD == FALSE){
    x.max <- max(dates.vector)
  }else{
    x.max <- parse_date_time(x = force.x.max.YY.MM.DD, orders = "ymd", tz = "Etc/GMT-5")
  }
  
  pts.span <- interval(start = floor_date(x.min, unit = label.with), end = ceiling_date(x.max, unit = label.with))
  
  num.ticks <- ceiling(pts.span / period(num = step.size, units = label.with)) + 1
  
  x.ticks <- get.tick.locs(num.ticks = num.ticks, label.with = label.with, step.size = step.size, pts.span = pts.span)
  
  x.labs <- get.x.labs(label.with = label.with, x.ticks = x.ticks, offset.labs = label.btwn.ticks)
  
  x.lab.locs <- get.lab.locs(label.with = label.with, x.ticks = x.ticks, offset.labs = label.btwn.ticks)
  
  # x.min <- min(x.ticks)
  # x.max <- max(x.ticks)
  
  return(list(sample.dates = dates.vector,
              tick.locations = x.ticks, 
              label.text = x.labs, 
              label.locations = x.lab.locs, 
              min = x.min, max = x.max))
}

# modified from pop:
set.up.stacked.years.plot <- function(dates.vector, year.range = NA, month.range = NA, rainbow.v = 1, rainbow.alpha = .3){
  # can specify as in year.range <- c(2000, 2017) and month.range <- c("2000-01-01", "2001-1-01")
  # if don't specify, max and min year and all months are used.
  
  if (is.na(year.range)){
    year.range <- c(min(year(dates.vector)), max(year(dates.vector)))
  }
  year.labs <- seq(from = year.range[1], to = year.range[2], by = 1)
  year.ticks <- c(year.labs, year.range[2] + 1) - .5
  
  if (is.na(month.range)){
    month.range <- c("2000-01-01", "2001-1-01")
  }
  month.range <- parse_date_time(x = month.range, orders = "ymd", tz = "Etc/GMT-5")
  month.ticks <- c(paste(2000, 1:12, 1, sep = "-"), "2001-1-1")
  month.ticks <- parse_date_time(x = month.ticks, orders = "ymd", tz = "Etc/GMT-5")
  month.lab.locs <- paste(2000, 1:12, 15, sep = "-")
  month.lab.locs <- parse_date_time(x = month.lab.locs, orders = "ymd", tz = "Etc/GMT-5")
  month.labs <- month(month.lab.locs, label = T, abbr = T)
  
  plot(x = month.range, y = year.range, type = "n", xlim = month.range, ylim = c(year.range[1] - .5, year.range[2] + .5), ann = F, axes = F)
  
  # month on x axis (adjust "line =" when change plot size)
  axis(side = 1, at = month.ticks, labels = F, line = -.75)
  axis(side = 1, at = month.lab.locs, labels = month.labs, tick = F, line = -1.25)
  
  # year on y axis (adjust "line =" when change plot size)
  axis(side = 2, at = year.ticks, labels = F, line = -1.65)
  axis(side = 2, at = year.labs, tick = F, las = 2, line = -2)
  
  # shading by month
  month.colors <- adjustcolor(col = rainbow(n = 13, v = rainbow.v), alpha.f = rainbow.alpha)
  month.colors <- c(month.colors[8:1],month.colors[12:9])
  for (m in 1:length(month.labs)){
    rect(xleft = month.ticks[m], xright = month.ticks[m + 1], ybottom = min(year.ticks), ytop = max(year.ticks), col = month.colors[m], border = NA)
  }
  return(year.labs)
}

add.dates.to.stacked.plot <- function(date.vector, year.labs, point.type = 23, point.col = adjustcolor(col = "grey", alpha.f = .4), line.col = "black", point.cex = 1){
  for (y in 1:length(year.labs)){
    index <- which(year(date.vector) == year.labs[y])
    if(length(index) > 0){
      x.days <- day(date.vector[index])
      x.months <- month(date.vector[index])
      x.vals <- paste(2000, x.months, x.days, sep = "-")
      x.vals <- parse_date_time(x = x.vals, orders = "ymd", tz = "Etc/GMT-5")
      y.vals <- rep.int(year.labs[y], times = length(index))
      if (length(point.col) > 1){
        pt.col <- point.col[index]
      }else{
        pt.col <- point.col
      }
      if (length(line.col) > 1){
        ln.col <- line.col[index]
      }else{
        ln.col <- line.col
      }
      if (length(point.cex) > 1){
        pt.cex <- point.cex[index]
      }else{
        pt.cex <- point.cex
      }
      if (length(point.type) > 1){
        pt.type <- point.type[index]
      }else{
        pt.type <- point.type
      }
      points(x = x.vals, y = y.vals, pch = pt.type, bg = pt.col, col = ln.col, cex = pt.cex)
    }
  }
  # "reminder: x-axis is in year 2000"
}

add.text.to.stacked.plot <- function(date.vector, text.vector, year.labs, point.type = 23, text.col = "black", text.cex = 1, text.srt = 0){
  for (y in 1:length(year.labs)){
    index <- which(year(date.vector) == year.labs[y])
    if(length(index) > 0){
      x.text <- text.vector[index]
      x.days <- day(date.vector[index])
      x.months <- month(date.vector[index])
      x.vals <- paste(2000, x.months, x.days, sep = "-")
      x.vals <- parse_date_time(x = x.vals, orders = "ymd", tz = "Etc/GMT-5")
      y.vals <- rep.int(year.labs[y], times = length(index))
      text(x = x.vals, y = y.vals, labels = x.text, cex = text.cex, col = text.col, srt = text.srt)
    }
  }
  # "reminder: x-axis is in year 2000"
}

add.legend.point <- function(date.vector = c("8-15-2020", "8-25-2020"), description = "label text", year.decimal = 0, point.type = 23, point.col = adjustcolor(col = "grey", alpha.f = .4), line.col = "black", point.cex = 1){
  # example date.vector = c("8-15-2020", "8-25-2020"), then y axis is at 2020, x axis for point is at 8/15 and description is at 8/25
  date.vector <- parse_date_time(x = date.vector, orders = "mdy", tz = "Etc/GMT-5")
  
  # draw point
  x.days <- day(date.vector[1])
  x.months <- month(date.vector[1])
  x.vals <- paste(2000, x.months, x.days, sep = "-")
  x.vals <- parse_date_time(x = x.vals, orders = "ymd", tz = "Etc/GMT-5")
  y.vals <- year(date.vector[1]) + year.decimal
  points(x = x.vals, y = y.vals, pch = point.type, col = line.col, bg = point.col, cex = point.cex, xpd = NA)
  
  # write text
  x.days <- day(date.vector[2])
  x.months <- month(date.vector[2])
  x.vals <- paste(2000, x.months, x.days, sep = "-")
  x.vals <- parse_date_time(x = x.vals, orders = "ymd", tz = "Etc/GMT-5")
  y.vals <- year(date.vector[2]) + year.decimal
  text(x = x.vals, y = y.vals, labels = description, adj = 0, xpd = NA)
  
  # reminder: x-axis is in year 2000
}


# ---- quick looks at assembly stats ----

hist(asms$Assembly.Size)
hist(asms$Num.Scaffolds)
hist(asms$Num.Genes)
hist(asms$Num.Bins)

# why are there 2 clear groups for these????
plot(asms$Num.Scaffolds ~ asms$Assembly.Size)
plot(asms$Num.Genes ~ asms$Assembly.Size)
plot(asms$Num.Bins ~ asms$Assembly.Size)
plot(asms$Num.Bins ~ asms$Num.Scaffolds)
plot(asms$Num.Bins ~ asms$Num.Genes)

asms$color <- "grey"
asms$color[which(asms$Sequencing.Plate == "tubes")] <- "red"
asms$color[which(asms$Sequencing.Plate == "plate.1")] <- "blue"
asms$color[which(asms$Sequencing.Plate == "plate.2")] <- "green"
asms$color[which(asms$Sequencing.Plate == "plate.3")] <- "orange"
asms$color[which(asms$Sequencing.Plate == "plate.4")] <- "purple"
asms$color[which(asms$Sequencing.Plate == "plate.5")] <- "pink"

plot(asms$Num.Scaffolds ~ asms$Assembly.Size, col = asms$color)
plot(asms$Num.Genes ~ asms$Assembly.Size, col = asms$color)
plot(asms$Num.Bins ~ asms$Assembly.Size, col = asms$color)
plot(asms$Num.Bins ~ asms$Num.Scaffolds, col = asms$color)
plot(asms$Num.Bins ~ asms$Num.Genes, col = asms$color)

# OK make one look decent and ask Neha
par(mfrow = c(2,2), mar = c(5,5,1,1), oma = c(.1,.1,3,.1))
plot(asms$Num.Scaffolds ~ asms$Assembly.Size, col = asms$color, ylab = "Number Scaffolds", xlab = "Genome Size")
plot(asms$Num.Genes ~ asms$Assembly.Size, col = asms$color, ylab = "Number Genes", xlab = "Genome Size")
plot(asms$Num.Bins ~ asms$Assembly.Size, col = asms$color, ylab = "Number Bins", xlab = "Genome Size")
plot(asms$Num.Bins ~ asms$Num.Scaffolds, col = asms$color, ylab = "Number Bins", xlab = "Number Scaffolds")
mtext(text = c("plate.1", "plate.2","plate.3","plate.4","plate.5","tubes"), side = 3, line = .75, 
      col = c("blue","green","orange","purple","pink","red"), outer = T, at = c(.2,.3,.4,.5,.6,.7), cex = 1.2)


# ---- quick look at overall bin taxonomies ----

colnames(bins)
for (c in 16:29){
  pdf(file = paste0("figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/Rank_Abund_",colnames(bins)[c],".pdf"), 
      width = 6, height = 5)
  par(mar = c(10,4,1,0))
  plot.bin.rank.abund(bins = bins, tax.col = c)
  dev.off()
  pdf(file = paste0("figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/Rank_Abund_",colnames(bins)[c],"-top_20.pdf"), 
      width = 6, height = 5 )
  par(mar = c(10,4,1,0))
  plot.bin.rank.abund(bins = bins, tax.col = c, nbar = 20)
  dev.off()
}

# ---- look at bins per metaG by stacked years ----
# BUBBLE DIAMETER REPRESENTS NUMBER BINS

sample.dates <- convert.colnames.to.dates(sample.names = asms$Sample.Name)
index <- order(sample.dates)
asms <- asms[index, ]

# text
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years.pdf",
    width = 40, height = 8)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = 21, point.col = "white", line.col = "white", point.cex = 3)
add.text.to.stacked.plot(date.vector = sample.dates, text.vector = asms$Num.Bins, year.labs = save.year.info, text.cex = .5, text.srt = 90)
dev.off()

hist(asms$Num.Bins)
point.col <- rep("white",nrow(asms))
point.pch <- rep(21,nrow(asms))
index <- which(asms$Num.Bins== 0)
point.col[index] <- "grey"
point.pch[index] <- 8

# text and greyed out zeros
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-grey_zeros-plate_outlines.pdf",
    width = 40, height = 8)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = 21, point.col = point.col, line.col = asms$color, point.cex = 3)
add.text.to.stacked.plot(date.vector = sample.dates, text.vector = asms$Num.Bins, year.labs = save.year.info, text.cex = .5, text.srt = 90)
mtext(text = c("plate.1", "plate.2","plate.3","plate.4","plate.5","tubes"), side = 3, line = -2, 
      col = c("blue","green","orange","purple","pink","red"), outer = T, at = c(.2,.3,.4,.5,.6,.7), cex = 1)
dev.off()

hist(asms$Num.Bins)
cex.vect <- rep(1,nrow(asms))
cex.vect <- cex.vect + asms$Num.Bins /25
hist(cex.vect)

# bubble and plate colors
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-bubble_num_bins-plate_colors.pdf",
    width = 40, height = 8)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = 21, point.col = adjustcolor(asms$color,.5), line.col = asms$color, point.cex = cex.vect)
mtext(text = c("plate.1", "plate.2","plate.3","plate.4","plate.5","tubes"), side = 3, line = -2, 
      col = c("blue","green","orange","purple","pink","red"), outer = T, at = c(.2,.3,.4,.5,.6,.7), cex = 1)
dev.off()

bin.size.legend <- c(0,10,25,40,80)
bin.cex.legend <- 1 + bin.size.legend /25

# bubbles and zeros grey
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-bubble_num_bins-zeros_grey.pdf",
    width = 40, height = 8)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = 21, point.col = adjustcolor(point.col,.5), line.col = point.col, point.cex = cex.vect)
rect(xleft = parse_date_time("2-1-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("9-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-20-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nBins", xpd = T)
add.legend.point(date.vector = c("3-18-2020","3-25-2020"), year.decimal = .65, description = bin.size.legend[1], point.type = 21, point.col = adjustcolor("grey",.5), line.col = "grey", point.cex = bin.cex.legend[1])
add.legend.point(date.vector = c("4-13-2020","4-20-2020"), year.decimal = .65, description = bin.size.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[2])
add.legend.point(date.vector = c("5-19-2020","5-29-2020"), year.decimal = .65, description = bin.size.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[3])
add.legend.point(date.vector = c("6-27-2020","7-9-2020"), year.decimal = .65, description = bin.size.legend[4], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[4])
add.legend.point(date.vector = c("8-10-2020","8-26-2020"), year.decimal = .65, description = bin.size.legend[5], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[5])
dev.off()

# bubbles and zeros grey - small plot
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-bubble_num_bins-zeros_grey-small_plot.pdf",
    width = 9, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = 21, point.col = adjustcolor(point.col,.5), line.col = point.col, point.cex = cex.vect)
rect(xleft = parse_date_time("2-1-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("9-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-20-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nBins", xpd = T)
add.legend.point(date.vector = c("3-18-2020","3-25-2020"), year.decimal = .65, description = bin.size.legend[1], point.type = 21, point.col = adjustcolor("grey",.5), line.col = "grey", point.cex = bin.cex.legend[1])
add.legend.point(date.vector = c("4-13-2020","4-20-2020"), year.decimal = .65, description = bin.size.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[2])
add.legend.point(date.vector = c("5-19-2020","5-29-2020"), year.decimal = .65, description = bin.size.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[3])
add.legend.point(date.vector = c("6-27-2020","7-9-2020"), year.decimal = .65, description = bin.size.legend[4], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[4])
add.legend.point(date.vector = c("8-10-2020","8-26-2020"), year.decimal = .65, description = bin.size.legend[5], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[5])
dev.off()

# bubbles, zeros grey, outline is plate color, small plot
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-bubble_num_bins-zeros_grey-plate_outline-small_plot.pdf",
    width = 9, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = 21, point.col = adjustcolor(point.col,.5), line.col = asms$color, point.cex = cex.vect)
rect(xleft = parse_date_time("2-1-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("9-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-20-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nBins", xpd = T)
add.legend.point(date.vector = c("3-18-2020","3-25-2020"), year.decimal = .65, description = bin.size.legend[1], point.type = 21, point.col = adjustcolor("grey",.5), line.col = "grey", point.cex = bin.cex.legend[1])
add.legend.point(date.vector = c("4-13-2020","4-20-2020"), year.decimal = .65, description = bin.size.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[2])
add.legend.point(date.vector = c("5-19-2020","5-29-2020"), year.decimal = .65, description = bin.size.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[3])
add.legend.point(date.vector = c("6-27-2020","7-9-2020"), year.decimal = .65, description = bin.size.legend[4], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[4])
add.legend.point(date.vector = c("8-10-2020","8-26-2020"), year.decimal = .65, description = bin.size.legend[5], point.type = 21, point.col = adjustcolor("white",.5), line.col = "white", point.cex = bin.cex.legend[5])
dev.off()






# bubbles and zeros asterisks - small plot
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-bubble_num_bins-zeros_asterisks-small_plot.pdf",
    width = 9, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = point.pch, 
                          point.col = adjustcolor("white",.5), line.col = "black", point.cex = cex.vect)
rect(xleft = parse_date_time("2-1-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("9-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-20-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nBins", xpd = T)
add.legend.point(date.vector = c("3-18-2020","3-25-2020"), year.decimal = .65, description = bin.size.legend[1], point.type = 8, point.col = adjustcolor("grey",.5), line.col = "black", point.cex = bin.cex.legend[1])
add.legend.point(date.vector = c("4-13-2020","4-20-2020"), year.decimal = .65, description = bin.size.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[2])
add.legend.point(date.vector = c("5-19-2020","5-29-2020"), year.decimal = .65, description = bin.size.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[3])
add.legend.point(date.vector = c("6-27-2020","7-9-2020"), year.decimal = .65, description = bin.size.legend[4], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[4])
add.legend.point(date.vector = c("8-10-2020","8-26-2020"), year.decimal = .65, description = bin.size.legend[5], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[5])
dev.off()

# bubbles and zeros asterisks - large plot
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-bubble_num_bins-zeros_asterisks.pdf",
    width = 40, height = 8)
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = point.pch, 
                          point.col = adjustcolor("white",.5), line.col = "black", point.cex = cex.vect)
rect(xleft = parse_date_time("2-1-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("9-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-20-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nBins", xpd = T)
add.legend.point(date.vector = c("3-18-2020","3-25-2020"), year.decimal = .65, description = bin.size.legend[1], point.type = 8, point.col = adjustcolor("grey",.5), line.col = "black", point.cex = bin.cex.legend[1])
add.legend.point(date.vector = c("4-13-2020","4-20-2020"), year.decimal = .65, description = bin.size.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[2])
add.legend.point(date.vector = c("5-19-2020","5-29-2020"), year.decimal = .65, description = bin.size.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[3])
add.legend.point(date.vector = c("6-27-2020","7-9-2020"), year.decimal = .65, description = bin.size.legend[4], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[4])
add.legend.point(date.vector = c("8-10-2020","8-26-2020"), year.decimal = .65, description = bin.size.legend[5], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[5])
dev.off()

# bubbles and zeros asterisks and plate outlines - small plot
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-bubble_num_bins-zeros_asterisks-outlines_plate-small_plot.pdf",
    width = 9, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = point.pch, 
                          point.col = adjustcolor("white",.5), line.col = asms$color, point.cex = cex.vect)
rect(xleft = parse_date_time("2-1-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("9-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-20-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nBins", xpd = T)
add.legend.point(date.vector = c("3-18-2020","3-25-2020"), year.decimal = .65, description = bin.size.legend[1], point.type = 8, point.col = adjustcolor("grey",.5), line.col = "black", point.cex = bin.cex.legend[1])
add.legend.point(date.vector = c("4-13-2020","4-20-2020"), year.decimal = .65, description = bin.size.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[2])
add.legend.point(date.vector = c("5-19-2020","5-29-2020"), year.decimal = .65, description = bin.size.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[3])
add.legend.point(date.vector = c("6-27-2020","7-9-2020"), year.decimal = .65, description = bin.size.legend[4], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[4])
add.legend.point(date.vector = c("8-10-2020","8-26-2020"), year.decimal = .65, description = bin.size.legend[5], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[5])
dev.off()

# bubbles and zeros asterisks and plate outlines - large plot
pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/bins_by_season-stacked_years-bubble_num_bins-zeros_grey.pdf",
    width = 40, height = 8)
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = point.pch, 
                          point.col = adjustcolor("white",.5), line.col = asms$color, point.cex = cex.vect)
rect(xleft = parse_date_time("2-1-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("9-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-20-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nBins", xpd = T)
add.legend.point(date.vector = c("3-18-2020","3-25-2020"), year.decimal = .65, description = bin.size.legend[1], point.type = 8, point.col = adjustcolor("grey",.5), line.col = "black", point.cex = bin.cex.legend[1])
add.legend.point(date.vector = c("4-13-2020","4-20-2020"), year.decimal = .65, description = bin.size.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[2])
add.legend.point(date.vector = c("5-19-2020","5-29-2020"), year.decimal = .65, description = bin.size.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[3])
add.legend.point(date.vector = c("6-27-2020","7-9-2020"), year.decimal = .65, description = bin.size.legend[4], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[4])
add.legend.point(date.vector = c("8-10-2020","8-26-2020"), year.decimal = .65, description = bin.size.legend[5], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bin.cex.legend[5])
dev.off()

# ---- look at genome size by stacked year ----
# BUBBLE AREA REPRESENTS GENOME SIZE

# OOH. OK someone answered my stack overflow :)
# cex scales the radius, not the area

colnames(asms)
options(scipen = 999)
hist(asms$Assembly.Size)
cex.vect <- rep(1,nrow(asms))
genome.scaled.down <- asms$Assembly.Size / (10^8) # call this the area
genome.cex <-  sqrt(genome.scaled.down / pi) # convert to "radius" to have dots scale with area
hist(genome.cex)

genome.pch <- rep(21,nrow(asms))
index <- which(asms$Assembly.Size == 0)
genome.pch[index] <- 8
genome.cex[index] <- 1

options(scipen = -1)
summary(asms$Assembly.Size)
bbp.legend <- c(7e8,1e9,2e09) 
bbp.cex <- bbp.legend / 10^8
bbp.cex <- sqrt(bbp.cex / pi)

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/genome_size-stacked_years-bubble-small_plot.pdf",
    width = 9, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = "black", point.cex = genome.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Assembly\nsize (bp)", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = bbp.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = bbp.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = bbp.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[3])
dev.off()

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/genome_size-stacked_years-bubble.pdf",
    width = 20, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = "black", point.cex = genome.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Assembly\nsize (bp)", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = bbp.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = bbp.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = bbp.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[3])
dev.off()

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/genome_size-stacked_years-bubble_plate-colors.pdf",
    width = 20, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = asms$color, point.cex = genome.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Assembly\nsize (bp)", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = bbp.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = bbp.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = bbp.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = bbp.cex[3])
dev.off()

# ---- look at number scaffolds by stacked year ----
# BUBBLE AREA SCALES WITH NUM SCAFFOLDS

colnames(asms)
hist(asms$Num.Scaffolds)
cex.vect <- rep(1,nrow(asms))
scaffold.cex <- asms$Num.Scaffolds / (10^5) # call this the bubble area
scaffold.cex <- sqrt(scaffold.cex/pi) # cex is scaling the radius
hist(scaffold.cex)

genome.pch <- rep(21,nrow(asms))
index <- which(asms$Num.Scaffolds == 0)
genome.pch[index] <- 8
scaffold.cex[index] <- 1

options(scipen = -1)
summary(asms$Num.Scaffolds)
scaf.legend <- c(5e5,10e5,3e6) 
scaf.cex <- sqrt( (scaf.legend / 10^5) / pi)

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/num_scaffolds-stacked_years-bubble-small_plot.pdf",
    width = 9, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaffold.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Number\nScaffolds", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = scaf.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = scaf.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = scaf.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[3])
dev.off()

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/num_scaffolds-stacked_years-bubble.pdf",
    width = 20, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaffold.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Number\nScaffolds", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = scaf.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = scaf.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = scaf.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[3])
dev.off()

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/num_scaffolds-stacked_years-bubble_plate-colors.pdf",
    width = 20, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = asms$color, point.cex = scaffold.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Number\nScaffolds", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = scaf.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = scaf.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = scaf.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = scaf.cex[3])
dev.off()

# ---- look at number genes by stacked year ----
# BUBBLE AREA SCALES WITH NUM Genes

colnames(asms)
hist(asms$Num.Genes)
genes.cex <- asms$Num.Scaffolds / (10^5) # call this the bubble area
genes.cex <- sqrt(genes.cex/pi) # cex is scaling the radius
hist(genes.cex)

genome.pch <- rep(21,nrow(asms))
index <- which(asms$Num.Scaffolds == 0)
genome.pch[index] <- 8
genes.cex[index] <- 1

options(scipen = -1)
summary(asms$Num.Genes)
genes.legend <- c(1e6,2e6,3e6) 
legend.cex <- sqrt( (genes.legend / 10^5) / pi)

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/num_genes-stacked_years-bubble-small_plot.pdf",
    width = 9, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = "black", point.cex = genes.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Number\nGenes", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = genes.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = genes.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = genes.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[3])
dev.off()

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/num_genes-stacked_years-bubble.pdf",
    width = 20, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = "black", point.cex = genes.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nGenes", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = genes.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = genes.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = genes.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[3])
dev.off()

pdf(file = "figures/2021-05-24_quick_looks_at_bins_and_assembly_stats/num_genes-stacked_years-bubble_plate-colors.pdf",
    width = 20, height = 6)
par(mar = c(3,3,2,0))
save.year.info <- set.up.stacked.years.plot(dates.vector = sample.dates, rainbow.v = .8, rainbow.alpha = .5)
add.dates.to.stacked.plot(date.vector = sample.dates, year.labs = save.year.info, point.type = genome.pch, 
                          point.col = adjustcolor("white",.5), line.col = asms$color, point.cex = genes.cex)
rect(xleft = parse_date_time("1-15-2000",orders = "mdy", tz = "Etc/GMT-5"), 
     xright = parse_date_time("10-15-2000",orders = "mdy", tz = "Etc/GMT-5"),
     ybottom = 2019.75, ytop = 2021.5, xpd = T, col = adjustcolor(rainbow(13,v = .8)[6],.5))
text(x = parse_date_time("2-25-2000",orders = "mdy", tz = "Etc/GMT-5"), y = 2020.65, labels = "Num\nGenes", xpd = T)
add.legend.point(date.vector = c("4-10-2020","4-20-2020"), year.decimal = .65, description = genes.legend[1], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[1])
add.legend.point(date.vector = c("6-10-2020","6-20-2020"), year.decimal = .65, description = genes.legend[2], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[2])
add.legend.point(date.vector = c("8-12-2020","8-24-2020"), year.decimal = .65, description = genes.legend[3], point.type = 21, point.col = adjustcolor("white",.5), line.col = "black", point.cex = legend.cex[3])
dev.off()
