# RRR
# need to request additional allocation
# look at exactly how much time we need

library(readxl)

tracker <- read_excel(path = "data/2022-07-07_tracking_mapping/tracking_jobs.xlsx")

plot.folder <- "figures/2022-08-08_how_long_did_mapping_sorting_take"

make.empty.list.structure <- function(ListNames){
  # the ListNames can be something like c("OTU", "kingdom","phylum","class","order","family/lineage","genus/clade","species/tribe")
  empty.list <- list(NULL)
  for (e in 1:length(ListNames)){
    empty.list[[e]] <- 0
    names(empty.list)[e] <- ListNames[e]
  }
  return(empty.list)
}


map.time <- tracker$`2.run.t`
map.hrs <- as.numeric(substr(x = map.time, start = 1, stop = 2))
map.min <- as.numeric(substr(x = map.time, start = 4, stop = 5))
map.sec <- as.numeric(substr(x = map.time, start = 7, stop = 8))
map.time <- (map.hrs * 60 * 60) + (map.min * 60) + map.sec # time is in seconds

map.N <- tracker$`2.N`
group <- tracker$Mapping.Group

map.SU.per.pair <- make.empty.list.structure(ListNames = unique(group))
for (g in unique(group)){
  index <- which(group == g)
  
  tot.pairs <- length(index)
  
  my.time <- map.time[index] / 60 # time in minutes
  my.nodes <- map.N[index]
  
  map.SU.per.pair[[g]] <- my.time * my.nodes * (1/tot.pairs)
  
}

boxplot(x = map.SU.per.pair) # y axis is node minutes per mapping pair


sort.time <- tracker$`3.run.t`
sort.hrs <- as.numeric(substr(x = sort.time, start = 1, stop = 2))
sort.min <- as.numeric(substr(x = sort.time, start = 4, stop = 5))
sort.sec <- as.numeric(substr(x = sort.time, start = 7, stop = 8))
sort.time <- (sort.hrs * 60 * 60) + (sort.min * 60) + sort.sec # time is in seconds

sort.N <- tracker$`3.N`
group <- tracker$Mapping.Group

sort.SU.per.pair <- make.empty.list.structure(ListNames = unique(group))
for (g in unique(group)){
  index <- which(group == g)
  
  tot.pairs <- length(index)
  
  my.time <- sort.time[index] / 60 # time in minutes
  my.nodes <- sort.N[index]
  
  sort.SU.per.pair[[g]] <- my.time * my.nodes * (1/tot.pairs)
  
}


boxplot(x = sort.SU.per.pair) # y axis is node minutes per mapping pair

# ----

summary(unlist(map.SU.per.pair))
sd(unlist(map.SU.per.pair), na.rm = T)

summary(unlist(sort.SU.per.pair))
sd(unlist(sort.SU.per.pair), na.rm = T)

pdf(file = file.path(plot.folder,"tacc_times.pdf"), width = 6.5, height = 2.5)

par(mfrow = c(1,2), mar = c(4,2,2,1), oma = c(.1,2.5,.1,.1))
boxplot(x = map.SU.per.pair, lty = 1, range = 0, axes = F, ann = F, col = "red3")
box()
axis(side = 1, at = 1:9, lwd = 0, lwd.ticks = 1, labels = F)
axis(side = 1, at = 1:9, lwd = 0, labels = c("test", "group 1", "group 2","group 3","group 4","group 5","group 6","group 7","group 8"), las = 2, line = -.25)
axis(side = 2, labels = F, lwd = 0, lwd.ticks = 1)
axis(side = 2, labels = T, lwd = 0, line = -.25, las = 2)
mtext(text = "Time per cross-mapping pair\n(node-minutes)", side = 2, line = .5, outer = T)
mtext(text = "Mapping", side = 3, line = .5)

boxplot(x = sort.SU.per.pair, lty = 1, range = 0, axes = F, ann = F, col = "red3")
box()
axis(side = 1, at = 1:9, lwd = 0, lwd.ticks = 1, labels = F)
axis(side = 1, at = 1:9, lwd = 0, labels = c("test", "group 1", "group 2","group 3","group 4","group 5","group 6","group 7","group 8"), las = 2, line = -.25)
axis(side = 2, labels = F, lwd = 0, lwd.ticks = 1)
axis(side = 2, labels = T, lwd = 0, line = -.25, las = 2)
mtext(text = "Sorting", side = 3, line = .5)

dev.off()

par(mfrow = c(1,1), mar = c(3,5,2,1))
boxplot(x = list("mapping" = unlist(map.SU.per.pair), "sorting" = unlist(sort.SU.per.pair)), lty = 1, range = 0, las = 2)

png(filename = file.path(plot.folder,"mapping_times.png"), width = 6.5, height = 2.5, units = "in", res = 300)
par(mfrow = c(1,1), mar = c(4,2,2,1), oma = c(.1,2.5,.1,.1))
boxplot(x = map.SU.per.pair, lty = 1, range = 0, axes = F, ann = F, col = "red3")
box()
axis(side = 1, at = 1:9, lwd = 0, lwd.ticks = 1, labels = F)
axis(side = 1, at = 1:9, lwd = 0, labels = c("test", "group 1", "group 2","group 3","group 4","group 5","group 6","group 7","group 8"), las = 2, line = -.25)
axis(side = 2, labels = F, lwd = 0, lwd.ticks = 1)
axis(side = 2, labels = T, lwd = 0, line = -.25, las = 2)
mtext(text = "Time per cross-mapping pair\n(node-minutes)", side = 2, line = .5, outer = T)
mtext(text = "Mapping", side = 3, line = .5)
dev.off()

png(filename = file.path(plot.folder,"sorting_times.png"), width = 6.5, height = 2.5, units = "in", res = 300)
par(mfrow = c(1,1), mar = c(4,2,2,1), oma = c(.1,2.5,.1,.1))
boxplot(x = sort.SU.per.pair, lty = 1, range = 0, axes = F, ann = F, col = "red3")
box()
axis(side = 1, at = 1:9, lwd = 0, lwd.ticks = 1, labels = F)
axis(side = 1, at = 1:9, lwd = 0, labels = c("test", "group 1", "group 2","group 3","group 4","group 5","group 6","group 7","group 8"), las = 2, line = -.25)
axis(side = 2, labels = F, lwd = 0, lwd.ticks = 1)
axis(side = 2, labels = T, lwd = 0, line = -.25, las = 2)
mtext(text = "Sorting", side = 3, line = .5)
dev.off()