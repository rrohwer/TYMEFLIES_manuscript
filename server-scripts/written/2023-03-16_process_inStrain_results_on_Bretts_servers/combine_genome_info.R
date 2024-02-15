# RRR
# read in the genome info tables from instrain
# combine them into a single file and save

# ---- set up ----

library(lubridate)
library(data.table)

user.input <- commandArgs(trailingOnly = TRUE)
folder.path <- user.input[1]    # ../../yggshare/current_projects/TYMEFLIES/tymeflies/runinstrain96/
sample.key <- user.input[2]     # drep_results_all_genomes_0.96.rds
season.key <- user.input[3]     # limony_seasons.rds
num.threads <- as.numeric(user.input[4])    # 40

# ---- functions ----

get.date.info <- function(sample.name, seasons, key){
  # get sample date and associated info for the genome_info file
  
  sample.name <- sub(pattern = "^.*output/*", replacement = "", x = sample.name)
  sample.name <- sub(pattern = ".IS_genome_info.tsv", replacement = "", x = sample.name)
  sample.date <- parse_date_time(x = substr(x = sample.name, start = 3, stop = 12), orders = "ymd", tz = "Etc/GMT-5")
  
  yr <- year(sample.date)
  seas <- seasons[seasons$Year == yr,-1] |>
    as.matrix() |>
    t() |>
    parse_date_time(orders = "ymd", tz = "Etc/GMT-5")
  seas <- data.frame("season" = colnames(seasons)[-1],"start" = seas)
  # there are no fall sample dates after the new year, so can take this shortcut for spring ice-on dates:
  seas <- rbind(data.frame("season" = "Ice.On", "start" = parse_date_time(x = paste(yr,1,1), orders = "ymd",tz = "Etc/GMT-5")), seas)
  
  index <- which(sample.date >= seas$start)
  my.season <- seas$season[max(index)]
  # my.season <- factor(x = my.season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")) # no need here b/c fwrite
  
  invasion <- "none"
  if(yr > 2009){
    invasion = "spiny"
  }
  if(yr > 2015){
    invasion = "zebra"
  }
  # invasion <- factor(x = invasion, levels = c("none","spiny","zebra")) # no need here b/c fwrite
  
  sample.year <- year(sample.date)
  sample.yday <- yday(sample.date)
  sample.date <- as.character(sample.date) # for faster fwrite
  
  return(list("name" = sample.name,"date" = sample.date, "year" = sample.year, "yday" = sample.yday,
              "season" = my.season, "invasion" = invasion))
}

# ---- get file paths ----
seasons <- readRDS(season.key)

sample.names <- readRDS(sample.key)
sample.names <- unique(sample.names$tymeflies.name)
my.files <- list.files(path = paste0(folder.path,"/",sample.names,".IS/","output/"), 
                       pattern = "genome_info", full.names = TRUE)

# ---- combine files ----

genome.info.list <- list()
for (f in my.files){
  one.genome <- fread(file = f, sep = "\t", header = T, nThread = num.threads)
  date.info <- get.date.info(sample.name = f, seasons = seasons)
  one.genome$sample <- date.info$name
  one.genome$date <- date.info$date
  one.genome$year <- date.info$year
  one.genome$yday <- date.info$yday
  one.genome$season <- date.info$season
  one.genome$invasion <- date.info$invasion
  
  genome.info.list <- c(genome.info.list, list(one.genome))
}

genome.info.table <- rbindlist(l = genome.info.list, use.names = TRUE) # **** really weird that cols not all in same order ****
rm(genome.info.list)

# ---- save output ----

fwrite(x = genome.info.table, file = "genome_info_combined.tsv.gz", sep = "\t", nThread = num.threads, compress = "gzip")

