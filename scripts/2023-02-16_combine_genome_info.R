# RRR
# read in the genome info tables from instrain
# combine them into a single file and save

library(lubridate)
library(limony)
data("seasons")

folder.path <- "data/2023-12-12_instrain_results_TIMEOUT"

my.files <- list.files(path = file.path(folder.path,"genome_info_files"), full.names = T)

get.date.info <- function(sample.name, seasons, key){
  sample.name <- sub(pattern = "^.*genome_info_files/", replacement = "", x = sample.name)
  sample.name <- sub(pattern = ".IS_genome_info.tsv", replacement = "", x = sample.name)
  
  sample.date <- parse_date_time(x = substr(x = sample.name, start = 3, stop = 12), orders = "ymd", tz = "Etc/GMT-5")
  
  yr <- year(sample.date)
  seas <- seasons[seasons$Year == yr,-1] |>
    as.matrix() |>
    t() |>
    parse_date_time(orders = "ymd", tz = "Etc/GMT-5")
  seas <- data.frame("season" = colnames(seasons)[-1],"start" = seas, "color" = unique(key$Color.Season))
  # there are no fall sample dates after the new year, so can take this shortcut for spring ice-on dates:
  seas <- rbind(data.frame("season" = "Ice.On",
                           "start" = parse_date_time(x = paste(yr,1,1), orders = "ymd",tz = "Etc/GMT-5"), 
                           "color" = "snow3"),
                seas)
  
  index <- which(sample.date >= seas$start)
  my.season <- seas$season[max(index)]
  my.season <- factor(x = my.season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))
  my.season.color <- factor(x = seas$color[max(index)], levels = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))
  
  invasion <- "none"
  invasion.color <- unique(key$Color.Invasion)[1]
  if(yr > 2009){
    invasion = "spiny"
    invasion.color = unique(key$Color.Invasion)[2]
  }
  if(yr > 2015){
    invasion = "zebra"
    invasion.color = unique(key$Color.Invasion)[3]
  }
  invasion <- factor(x = invasion, levels = c("none","spiny","zebra"))
  invasion.color <- factor(x = invasion.color, levels = c("steelblue","orange2","red3"))
  
  return(list("name" = sample.name,"date" = sample.date, 
              "season" = my.season, "color.season" = my.season.color, 
              "invasion" = invasion, "color.invasion" = invasion.color))
}

genome.info <- read.table(file = my.files[1], sep = "\t", header = T)
date.info <- get.date.info(sample.name = my.files[1], seasons = seasons, key = key)
genome.info$sample <- date.info$name
genome.info$date <- date.info$date
genome.info$season <- date.info$season
genome.info$invasion <- date.info$invasion
genome.info$color.season <- date.info$color.season
genome.info$color.invasion <- date.info$color.invasion

for (f in 2:length(my.files)){
  one.genome <- read.table(file = my.files[f], sep = "\t", header = T)
  date.info <- get.date.info(sample.name = my.files[f], seasons = seasons, key = key)
  one.genome$sample <- date.info$name
  one.genome$date <- date.info$date
  one.genome$season <- date.info$season
  one.genome$invasion <- date.info$invasion
  one.genome$color.season <- date.info$color.season
  one.genome$color.invasion <- date.info$color.invasion
  
  genome.info <- rbind(genome.info, one.genome)
}

saveRDS(object = genome.info, file = file.path(folder.path,"genome_info_combined.rds"))
