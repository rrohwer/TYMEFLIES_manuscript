# RRR
library(viridis)
library(gt)
library(scales)
library(lubridate)
library(patchwork)
library(ggplot2)

bins <- readRDS("data/2022-09-23_bin_stats/all_bins.rds")
colnames(bins)
head(bins)

plot.folder <- "figures/2022-09-25_explore_the_bins/"

# ---- choose taxa ----

tax.lvl <- "order"
t.lvl <- which(colnames(bins) == tax.lvl)
taxon <- "o__Nanopelagicales"

tax.lvl <- "class"
t.lvl <- which(colnames(bins) == tax.lvl)
taxon <- "c__Cyanobacteriia"

tax.lvl <- "species"
t.lvl <- which(colnames(bins) == tax.lvl)
taxon <- "s__Dolichospermum flosaquae"

tax.lvl <- "family"
t.lvl <- which(colnames(bins) == tax.lvl)
taxon <- "f__Burkholderiaceae"

unique(bins$domain)
# no archaea in GTDB?
unique(bins$phylum) # can choose any of these
tax.lvl <- "phylum"
t.lvl <- which(colnames(bins) == tax.lvl)
taxon <- "p__Bacteroidota"
taxon <- "p__Actinobacteriota"
taxon <- "p__Proteobacteria"
taxon <- "p__Acidobacteriota"
taxon <- "p__Verrucomicrobiota"
taxon <- "p__Planctomycetota"
taxon <- "p__Chloroflexota"

index <- (which(bins$phylum == "p__4484-113" |
                  bins$phylum == "p__" |
                  bins$phylum == "p__J088" |
                  bins$phylum == "p__JAGOBX01" |
                  bins$phylum == "p__UBP6" ))
                  # bins$phylum == "p__SAR324" ))
my.org <- bins[index, ]
taxon <- "weird ones"

# ---- select org ----

index <- which(bins[,tax.lvl] == taxon)

my.org <- bins[index, ]

# or don't subset
my.org <- bins
taxon <- "All Bins"
t.lvl <- 0
tax.lvl <- "all"

# ---- plot name ----

plot.name <- paste(t.lvl,tax.lvl,taxon, sep = "_")

# ---- how many bins are there ----

hist(my.org$Completeness)
hist(my.org$Contamination)

num.bins <- matrix(data = NA, nrow = 5, ncol = 11)
completeness <- seq.int(from = 0, to = 100, by = 10)
redundancy <- c(1,2.5,5,10,100)
colnames(num.bins) <- paste("over",completeness, sep = "_")
row.names(num.bins) <- paste("under", redundancy, sep = "_")
for (c in 1:length(completeness)){
  for (r in 1:length(redundancy)){
    num.bins[r,c] <- sum(my.org$Completeness >= completeness[c] & my.org$Contamination <= redundancy[r])
  }
}
num.bins

my.tab <- gt(data = as.data.frame(num.bins[5:1, ]), rownames_to_stub = T, caption = taxon)
# data_color(data = my.tab, columns = 1:12, colors = col_numeric(palette = viridis(100,direction = -1), domain = 0:max(num.bins)))
my.tab <- tab_style(data = my.tab, style = cell_fill(color = "orange1"), locations = cells_body(columns = "over_50", rows = "under_10"))
my.tab <- tab_style(data = my.tab, style = cell_fill(color = "red"), locations = cells_body(columns = "over_90", rows = "under_5"))
my.tab

# ---- save table ----

gtsave(data = my.tab, filename = file.path(plot.folder,paste0(plot.name,"-table.html")))

# ---- look at other levels ----

MQ.completeness <- 50 # from https://www.nature.com/articles/nbt.3893/tables/1 but ignoring RNA requirements
MQ.contamination <- 10
HQ.completeness <- 90
HQ.contamination <- 5


for (t in (t.lvl):ncol(my.org)){
  daughter.names <- unique(my.org[ ,t])
  temp <- my.org
  temp$LQ <- my.org$Completeness < MQ.completeness | my.org$Contamination >= MQ.contamination
  temp$HQ <- my.org$Completeness >= HQ.completeness & my.org$Contamination < HQ.contamination
  temp$MQ <- my.org$Completeness >= MQ.completeness & my.org$Contamination < MQ.contamination & (!temp$HQ)
  cat(all(rowSums(x = temp[ ,c("LQ","MQ","HQ")]) == 1))
  my.abunds <- aggregate(x = temp[ ,c("LQ","MQ","HQ")], by = list(temp[ ,t]), FUN = sum)
  row.names(my.abunds) <- my.abunds[ ,1]
  my.abunds <- my.abunds[ ,-1]
  my.abunds <- as.matrix(my.abunds)
  tot.abunds <-rowSums(my.abunds)
  my.abunds <- my.abunds[order(tot.abunds, decreasing = F), ]
  my.abunds <- t(my.abunds)
  tot.abunds <- colSums(my.abunds)
  
  
  if (nrow(my.abunds) == 1){
    png(filename = file.path(plot.folder,paste0(plot.name,"-",t,".png")), width = 6.5, height = 6.5, units = "in", res = 300)
    par(mar = c(3,4,6,1))
    bar.spots <- barplot(height = my.abunds[1, ], col = c("grey","orange1","red"), axes = F)
    axis(side = 2, las = 2)
    text(x = bar.spots, y = my.abunds[1, ] + max(tot.abunds) / 25, labels = my.abunds[1, ], xpd = NA)
    mtext(text = paste("Number of", taxon, "bins"), cex = 1.5, line = 3)
    dev.off()
    next
  }
  
  png(filename = file.path(plot.folder,paste0(plot.name,"-",t,".png")), width = 6.5, height = 6.5, units = "in", res = 300)
  par(mar = c(3,15,3,7))
  bar.spots <- barplot(height = my.abunds, names.arg = colnames(my.abunds), las = 1, horiz = T, border = T, col = c("grey","orange1","red"))
  text(y = bar.spots, x = tot.abunds + max(tot.abunds)/30, labels = paste(my.abunds[1, ], my.abunds[2, ], my.abunds[3, ], sep = ", "), xpd = NA, adj = 0, cex = .7)
  mtext(text = paste(colnames(my.org)[t], "under",taxon), cex = 1.5)
  mtext(text = "(Number of LQ, MQ, HQ bins)", cex = .7, side = 1, line = -5, at = max(tot.abunds) - max(tot.abunds)/20)
  dev.off()
}

# Notes:
# Ca. Planktophila dulcis 	acI-A1		
# Ca. Planktophila sulfonica 	acI-A1
# Ca. Planktophila versatilis 	acI-A1
# Ca. Planktophila lacus 	acI-A1
# Ca. Planktophila rubra 	acI-A1
# Ca. Planktophila limnetica	acI-A2		
# Ca. Planktophila aquatilis 	acI-A4
# Ca. Planktophila vernalis 	acI-A7

# Ca. Planktophila limnetica 	Phila

# Ca. Nanopelagicus limnes 	acI-B1
# Ca. Nanopelagicus hibericus 	acI-B1
# Ca. Nanopelagicus abundans 	acI-B1

# is it possible the aphanizominon is being called dolichospermum? https://forum.gtdb.ecogenomic.org/t/about-the-bacterial-taxonomy-category/17
# yes you can see here it is called s_Dolichospermum flosaquae https://gtdb.ecogenomic.org/searches?s=al&q=aphanizomenon%20flos-aquae






# ---- and what is acI spread over time like? ----

fill.under.lines <- function(X, Y, YAxisMin, Color, xpd = F){
  poly.x <- c(min(X), X, max(X))
  poly.y <- c(YAxisMin, Y, YAxisMin )
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

sample.dates <- substr(x = my.org$Name, start = 3, stop = 12)
sample.dates <- parse_date_time(sample.dates, "ymd")

temp <- my.org
temp$LQ <- my.org$Completeness < MQ.completeness | my.org$Contamination >= MQ.contamination
temp$HQ <- my.org$Completeness >= HQ.completeness & my.org$Contamination < HQ.contamination
temp$MQ <- my.org$Completeness >= MQ.completeness & my.org$Contamination < MQ.contamination & (!temp$HQ)
cat(all(rowSums(x = temp[ ,c("LQ","MQ","HQ")]) == 1))
my.abunds <- aggregate(x = temp[ ,c("LQ","MQ","HQ")], by = list(sample.dates), FUN = sum)
colnames(my.abunds)[1] <- "date"
str(my.abunds)

par(mar = c(2,3,2,1))
plot(x = my.abunds$date, y = my.abunds$LQ, type = "n")
fill.under.lines(X = my.abunds$date, Y = my.abunds$LQ, YAxisMin = 0, Color = adjustcolor("grey",.3))
fill.under.lines(X = my.abunds$date, Y = my.abunds$MQ, YAxisMin = 0, Color = adjustcolor("orange1",.3))
fill.under.lines(X = my.abunds$date, Y = my.abunds$HQ, YAxisMin = 0, Color = adjustcolor("red",.3))
lines(x = my.abunds$date, y = my.abunds$LQ, col = "grey")
lines(x = my.abunds$date, y = my.abunds$MQ, col = "orange1")
lines(x = my.abunds$date, y = my.abunds$HQ, col = "red")

# cut out outlier
par(mar = c(2,3,2,1))
plot(x = my.abunds$date, y = my.abunds$LQ, type = "n", ylim = c(0,50))
fill.under.lines(X = my.abunds$date, Y = my.abunds$LQ, YAxisMin = 0, Color = adjustcolor("grey",.3))
fill.under.lines(X = my.abunds$date, Y = my.abunds$MQ, YAxisMin = 0, Color = adjustcolor("orange1",.3))
fill.under.lines(X = my.abunds$date, Y = my.abunds$HQ, YAxisMin = 0, Color = adjustcolor("red",.3))

lines(x = my.abunds$date, y = my.abunds$LQ, col = "grey")
lines(x = my.abunds$date, y = my.abunds$MQ, col = "orange1")
lines(x = my.abunds$date, y = my.abunds$HQ, col = "red")

points(x = my.abunds$date, y = my.abunds$LQ, col = "grey", cex = .5, pch = 19)
points(x = my.abunds$date, y = my.abunds$MQ, col = "orange1", cex = .5, pch = 19)
points(x = my.abunds$date, y = my.abunds$HQ, col = "red", cex = .5, pch = 19)

# interesting that the number of LQ bins mirrors the number of HQ bins, instead of being inverse

par(mfrow = c(10,2), mar = c(1,1,1,1), oma = c(1,1,2,0))
for (y in 2000:2019){
  my.year <- my.abunds[year(my.abunds$date) == y, ]
  plot(x = parse_date_time(c(paste(y,1,1),paste(y,12,31)),"ymd"), y = c(0,50), type = "n")
  fill.under.lines(X = my.year$date, Y = my.year$LQ, YAxisMin = 0, Color = adjustcolor("grey",.3))
  fill.under.lines(X = my.year$date, Y = my.year$MQ, YAxisMin = 0, Color = adjustcolor("orange1",.3))
  fill.under.lines(X = my.year$date, Y = my.year$HQ, YAxisMin = 0, Color = adjustcolor("red",.3))
  mtext(text = y, side = 2, line = -4, las = 2)
}
mtext(outer = T, text = taxon)

