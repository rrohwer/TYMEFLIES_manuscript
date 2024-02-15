# RRR 
# plot  sample dates over time for an overview of when the time-series data is from.
# scale or color sample dates by their attributes
# add lines to the plot, for example the secchi depths

# the input is a character-formatted date vector. Get it using:
#   import_metadata_excel_files.R to import the sample data from the excel sheets
#   use combine.date.columns() to format a character date vector "MM-DD-YYYY" of the samples you want
#   export date vector as .csv file with "row.names = F"

# ---- libraries ----

library("lubridate")

# ---- functions ----

convert.colnames.to.dates <- function(sample.names){
  sample.names <- substr(x = sample.names, start = 1, stop = 9)
  sample.names <- parse_date_time(x = sample.names, orders = "dmy", tz = "Etc/GMT-5")
  return(sample.names)
}

set.up.plot <- function(dates.vector = NA, year.range = NA, month.range = NA, shade.by.month = FALSE, year.lines = FALSE, month.lines = FALSE,
                        x.ax.lab.line = -1.25, x.ax.tck = -.025, x.ax.las = 1, y.ax.lab.line = -2, y.ax.tck = -.025, 
                        y.ax.btwn.lines = TRUE, x.ax.btwn.lines = TRUE, y.ax.tick.freq = 1, x.ax.tick.freq = 1, x.ax.single.letters = FALSE, y.ax.cex = 1, x.ax.cex = 1){
  # can specify as in year.range <- c(2000, 2017) and month.range <- c("2000-01-01", "2001-1-01")
  # if don't specify, max and min year and all months are used.
  
  if (is.na(year.range[1])){
    year.range <- c(min(year(dates.vector)), max(year(dates.vector)))
  }
  year.labs <- seq(from = year.range[1], to = year.range[2], by = 1)
  year.ticks <- c(year.labs, year.range[2] + 1) - .5
  year.labs <- year.labs[seq(from = 1, to = length(year.labs), by = y.ax.tick.freq)]
  year.ticks <- year.ticks[seq(from = 1, to = length(year.ticks), by = y.ax.tick.freq)]
  if(!y.ax.btwn.lines){
    year.ticks <- year.labs
  }
  
  if (is.na(month.range[1])){
    month.range <- c("2000-01-01", "2001-1-01")
  }
  month.range <- parse_date_time(x = month.range, orders = "ymd", tz = "Etc/GMT-5")
  month.ticks <- c(paste(2000, 1:12, 1, sep = "-"), "2001-1-1")
  month.ticks <- parse_date_time(x = month.ticks, orders = "ymd", tz = "Etc/GMT-5")
  month.ticks <- month.ticks[seq(from = 1, to = length(month.ticks), by = x.ax.tick.freq)]
  month.lab.locs <- paste(2000, 1:12, 15, sep = "-")
  month.lab.locs <- parse_date_time(x = month.lab.locs, orders = "ymd", tz = "Etc/GMT-5")
  month.lab.locs <- month.lab.locs[seq(from = 1, to = length(month.lab.locs), by = x.ax.tick.freq)]
  if(!x.ax.btwn.lines){
    month.ticks <- month.ticks[-length(month.ticks)]
    month.lab.locs <- month.ticks
  }
  month.labs <- lubridate::month(month.lab.locs, label = T, abbr = T)
  if (x.ax.single.letters){
    month.labs <- substr(x = month.labs, start = 1, stop = 1)
  }
  
  plot(x = month.range, y = year.range, type = "n", xlim = month.range, ylim = c(year.range[1] - .5, year.range[2] + .5), ann = F, axes = F, xaxs = "i", yaxs = "i")
  
  box()
  
  # month on x axis (adjust "line =" when change plot size)
  axis(side = 1, at = month.ticks, labels = F, line = 0, tck = x.ax.tck, lwd = 0, lwd.ticks = 1)
  axis(side = 1, at = month.lab.locs, labels = month.labs, lwd = 0, line = x.ax.lab.line, las = x.ax.las, cex.axis = x.ax.cex)
  
  # year on y axis (adjust "line =" when change plot size)
  axis(side = 2, at = year.ticks, labels = F, line = 0, tck = y.ax.tck, lwd = 0, lwd.ticks = 1)
  axis(side = 2, at = year.labs, lwd = 0, las = 2, line = y.ax.lab.line, cex.axis = y.ax.cex)
  
  # shading by month
  if (shade.by.month){
    month.colors <- adjustcolor(col = rainbow(n = 13), alpha.f = .3)
    month.colors <- c(month.colors[8:1],month.colors[12:9])
    for (m in 1:length(month.labs)){
      rect(xleft = month.ticks[m], xright = month.ticks[m + 1], ybottom = min(year.ticks), ytop = max(year.ticks), col = month.colors[m], border = NA)
    }
  }
  
  # year or month gridlines
  if (year.lines){
    abline(h = year.ticks, col = adjustcolor("black",alpha.f = .3), xpd = F)
  }
  if (month.lines){
    abline(v = month.ticks, col = adjustcolor("black", alpha.f = .3), xpd = F)
  }
  
  return(year.labs)
}

shade.by.season <- function(season.starts, season.ends, season.color, axis.year = 2000){
  no.NAs <- !is.na(season.starts) & !is.na(season.ends)
  season.starts <- season.starts[no.NAs]
  season.ends <- season.ends[no.NAs]
  
  year.vect <- year(season.starts)
  
  start.vect <- paste(axis.year, month(season.starts), day(season.starts))
  end.vect <-  paste(axis.year, month(season.ends), day(season.ends))
  
  start.vect <- parse_date_time(x = start.vect, orders = "ymd", tz = "Etc/GMT-5")
  end.vect <- parse_date_time(x = end.vect, orders = "ymd", tz = "Etc/GMT-5")
  
  for (y in 1:length(year.vect)){
    rect(xleft = start.vect[y], ybottom = year.vect[y] - .5, xright = end.vect[y], ytop = year.vect[y] + .5, 
         col = season.color, border = NA)
    
  }
}

get.year.scaling.and.offset <- function(my.max, my.min){
  my.range = my.max - my.min
  year.range = 1
  my.scaling <- year.range / my.range
  my.offset <- .5 - my.max * my.scaling
  return(list("scaling" = my.scaling, "offset" = my.offset))
}

get.new.x.y <- function(x, y, axis.year = 2000, year.offset = 0, year.scaling = 1){
  new.x <- paste(axis.year, month(x), day(x))
  new.x <- parse_date_time(x = new.x, orders = "ymd", tz = "Etc/GMT-5")
  data.year <- year(x)
  new.y <- y * year.scaling
  new.y <- data.year + year.offset + new.y
  return(list("x" = new.x, "y" = new.y, "year" = data.year))
}

get.scaled.point.cex <- function(min.cex, max.cex, my.values){
  unscaled.min <- min(my.values)
  unscaled.max <- max(my.values)
  
  # y = mx + b
  m = (max.cex - min.cex) / (unscaled.max - unscaled.min)
  b = min.cex - m * unscaled.min
  
  scaled.cex <- my.values * m + b
  return(scaled.cex)
}

add.dates <- function(date.vector, year.labs, point.type = 23, point.col = adjustcolor(col = "grey", alpha.f = .4), line.col = "black", point.cex = 1, offset = 0){
  # offset = .5 to be at top of year section or -.5 to be at bottom of year section
  for (y in 1:length(year.labs)){
    index <- which(year(date.vector) == year.labs[y])
    if(length(index) > 0){
      x.days <- day(date.vector[index])
      x.months <- month(date.vector[index])
      x.vals <- paste(2000, x.months, x.days, sep = "-")
      x.vals <- parse_date_time(x = x.vals, orders = "ymd", tz = "Etc/GMT-5")
      y.vals <- rep.int(year.labs[y], times = length(index))
      y.vals <- y.vals + offset
      points(x = x.vals, y = y.vals, pch = point.type, bg = point.col, col = line.col, cex = point.cex)
    }
  }
  # "reminder: x-axis is in year 2000"
}

add.text <- function(date.vector, year.labs, text.vector, text.col = "black", text.cex = 1, text.srt = 0, offset = 0){
  # offset = .5 to be at top of year section or -.5 to be at bottom of year section
  for (y in 1:length(year.labs)){
    index <- which(year(date.vector) == year.labs[y])
    if(length(index) > 0){
      x.days <- day(date.vector[index])
      x.months <- month(date.vector[index])
      x.vals <- paste(2000, x.months, x.days, sep = "-")
      x.vals <- parse_date_time(x = x.vals, orders = "ymd", tz = "Etc/GMT-5")
      y.vals <- rep.int(year.labs[y], times = length(index))
      y.vals <- y.vals + offset
      text(x = x.vals, y = y.vals, labels = text.vector[index], col = text.col, cex = text.cex, srt = text.srt)
    }
  }
  # "reminder: x-axis is in year 2000"
}

add.legend.point <- function(date.vector = c("8-15-2020", "8-25-2020"), description = "label text", point.type = 23, point.col = adjustcolor(col = "grey", alpha.f = .4), line.col = "black", point.cex = 1, point.vertical.scoot = 0, text.vertical.scoot = 0, text.col = "black", text.adj = 0, text.cex = 1)     {
  # example date.vector = c("8-15-2020", "8-25-2020"), then y axis is at 2020, x axis for point is at 8/15 and description is at 8/25
  date.vector <- parse_date_time(x = date.vector, orders = "mdy", tz = "Etc/GMT-5")
  
  # draw point
  x.days <- day(date.vector[1])
  x.months <- month(date.vector[1])
  x.vals <- paste(2000, x.months, x.days, sep = "-")
  x.vals <- parse_date_time(x = x.vals, orders = "ymd", tz = "Etc/GMT-5")
  y.vals <- year(date.vector[1]) + point.vertical.scoot
  points(x = x.vals, y = y.vals, pch = point.type, col = line.col, bg = point.col, cex = point.cex, xpd = NA)
  
  # write text
  x.days <- day(date.vector[2])
  x.months <- month(date.vector[2])
  x.vals <- paste(2000, x.months, x.days, sep = "-")
  x.vals <- parse_date_time(x = x.vals, orders = "ymd", tz = "Etc/GMT-5")
  y.vals <- year(date.vector[2]) + text.vertical.scoot
  text(x = x.vals, y = y.vals, labels = description, adj = text.adj, xpd = NA, col = text.col, cex = text.cex)
  
  # reminder: x-axis is in year 2000
}

add.side.legend.points <- function(year = c(2019), month.day = "1-15", point.type = 23, point.col = adjustcolor(col = "grey", alpha.f = .4), line.col = "black", point.cex = 1){
  
  x.loc <- parse_date_time(x = paste(2001, month.day), orders = "ymd", tz = "Etc/GMT-5")
  if(length(x.loc) < length(year)){
    x.loc <- rep(x.loc, length(year))
  }
  points(x = x.loc, y = year, pch = point.type, col = line.col, bg = point.col, cex = point.cex, xpd = NA)
  
  # reminder: x-axis is in year 2000
}

add.side.legend.text <- function(year = c(2019), month.day = "1-15", description = "label text", text.col = "black", text.adj = 0, text.cex = 1){
 
  x.loc <- parse_date_time(x = paste(2001, month.day), orders = "ymd", tz = "Etc/GMT-5")
  if(length(x.loc) < length(year)){
    x.loc <- rep(x.loc, length(year))
  }
  text(x = x.loc, y = year, labels = description, adj = text.adj, xpd = NA, col = text.col, cex = text.cex)
  
  # reminder: x-axis is in year 2000
}

# ---- 



