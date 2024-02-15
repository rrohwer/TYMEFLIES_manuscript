# RRR
library(lubridate)
library(readxl)

bigtyme <- read_excel("data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx")
yearly <- readRDS("pop/data/environmental_data/Robin-Refined/seasons/8-yearly_season_starts.rds")

split.by.year <- "data/2022-05-24_metabat_by_year/samplenames_by_year.rds"
split.by.subsample <- "data/2022-05-24_metabat_by_year/samplenames_by_year_subsampled.rds"
split.by.season <- "data/2022-05-24_metabat_by_year/samplenames_by_year_and_season.rds"

# ---- functions ----

make.empty.list.structure <- function(ListNames){
  # the ListNames can be something like c("OTU", "kingdom","phylum","class","order","family/lineage","genus/clade","species/tribe")
  empty.list <- list(NULL)
  for (e in 1:length(ListNames)){
    empty.list[[e]] <- 0
    names(empty.list)[e] <- ListNames[e]
  }
  return(empty.list)
}

# ---- go ----

bigtyme <- bigtyme[!is.na(bigtyme$taxon_oid), ]
bigtyme$sample.names <- paste0(bigtyme$TYMEFLIES.name, "_",bigtyme$taxon_oid)
bigtyme$sample.dates <- parse_date_time(x = paste(bigtyme$Year,bigtyme$Month,bigtyme$Day), orders = "ymd", tz = "Etc/GMT-5")
bigtyme <- bigtyme[order(bigtyme$sample.dates), ]

my.years <- sort(unique(bigtyme$Year))
samples.by.year <- make.empty.list.structure(ListNames = my.years)
num.per.year <- NULL
for (y in 1:length(my.years)){
  yr <- my.years[y]
  i <- which(bigtyme$Year == yr)
  num.per.year[y] <- length(i)
  names(num.per.year)[y] <- yr
  cat(yr,length(i),"\n")
  samples.by.year[[y]] <- bigtyme$sample.names[i]
}

est.hrs <- num.per.year * num.per.year * 40 /60
est.days <- est.hrs / 24
cbind(est.hrs, est.days,num.per.year)
colSums(cbind(est.hrs, est.days))
145 /30 # 4 months

saveRDS(object = samples.by.year, file = split.by.year)

hist(est.days, breaks = 20)
hist(num.per.year)

# ---- split up the long years by sub-sampling ----

index <- which(num.per.year > 20)
num.per.year[index]

samples.by.subsample <- make.empty.list.structure(ListNames = c(2000:2011,paste0(rep(2012:2018, each = 2),c("A","B")),2019))
num.per.subsample = NULL
for (y in 1:length(samples.by.subsample)){
  yr <- as.numeric(substr(x = names(samples.by.subsample), start = 1, stop = 4)[y])
  i <- which(bigtyme$Year == yr)
  do.ss <- length(i) > 20
  if (do.ss){
    ss <- substr(x = names(samples.by.subsample), start = 5, stop = 5)[y]
    if (ss == "A"){ 
      start.pos <- 1
    }else if (ss == "B"){
      start.pos <- 2
    }
    i <- i[seq.int(from = start.pos, to = length(i), by = 2)]
  }
  num.per.subsample[y] <- length(i)
  names(num.per.subsample)[y] <- yr
  cat(yr,length(i),"\n")
  samples.by.subsample[[y]] <- bigtyme$sample.names[i]
}

est.hrs <- num.per.subsample * num.per.subsample * 13 /60
est.days <- est.hrs / 24
cbind(est.hrs, est.days,num.per.subsample)
colSums(cbind(est.hrs, est.days))
84 /30 # 3 months

saveRDS(object = samples.by.subsample, file = split.by.subsample)

hist(est.days)
hist(num.per.subsample, breaks = 20)

# ---- split the long years by season ----

# I think this makes less sense... or maybe more sense, I dunno!

index <- which(num.per.year > 20)
num.per.year[index]

samples.by.season <- make.empty.list.structure(ListNames = c(2000:2011,paste0(rep(2012:2018, each = 2),c("A","B")),2019))
num.per.season = NULL
for (y in 1:length(samples.by.season)){
  yr <- as.numeric(substr(x = names(samples.by.season), start = 1, stop = 4)[y])
  seas <- yearly[yearly$Year == yr, ]
  i <- which(bigtyme$Year == yr)
  my.dates <- bigtyme$sample.dates[i]
  do.ss <- length(i) > 20
  if (do.ss){
    ss <- substr(x = names(samples.by.season), start = 5, stop = 5)[y]
    if (ss == "A"){ 
      i <- i[which(my.dates < seas$`Late Summer`)]
      bigtyme$sample.dates[i]
    }else if (ss == "B"){
      i <- i[which(my.dates >= seas$`Late Summer`)]
      bigtyme$sample.dates[i]
    }
  }
  num.per.season[y] <- length(i)
  names(num.per.season)[y] <- yr
  cat(yr,length(i),"\n")
  samples.by.season[[y]] <- bigtyme$sample.names[i]
}

est.hrs <- num.per.season * num.per.season * 13 /60
est.days <- est.hrs / 24
cbind(est.hrs, est.days,num.per.season)
colSums(cbind(est.hrs, est.days))
89 /30 # 3 months

saveRDS(object = samples.by.season, file = split.by.season)

hist(est.days)

# this takes longer b/c late summer + fall is more samples than half of the year. 
# could do fall to early summer, but gets complex splitting years like that.
# also, we need *differential* abundance signals, so subsampling full year might be better
# but if it's only IN that part of the year, subsampling seasons might be better...
# halp...




# ---- split into groups of 50ish ----

# Did this manually. Aim for 50 in a group, but also don't group 
