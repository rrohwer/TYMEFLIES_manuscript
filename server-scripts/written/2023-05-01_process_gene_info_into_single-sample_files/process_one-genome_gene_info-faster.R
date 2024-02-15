# RRR
# read in the gene info tables from instrain
# combine them into a single file and save
# edited the genome info script to go faster with data.table fread and rbindlist

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
  # get sample date and associated info for the gene_info file
  
  sample.name <- sub(pattern = "^.*output/*", replacement = "", x = sample.name)
  sample.name <- sub(pattern = ".IS_gene_info.$", replacement = "", x = sample.name)
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
                       pattern = "gene_info", full.names = TRUE) # some are .tsv and some are .tsv.gz ???

# ---- combine files ----

gene.info.list <- list()
for (f in my.files){
  cat("processing ",f,"\n")
  one.gene <- fread(file = f, sep = "\t", header = T, nThread = num.threads)
  date.info <- get.date.info(sample.name = f, seasons = seasons)
  one.gene$sample <- date.info$name
  one.gene$date <- date.info$date
  one.gene$year <- date.info$year
  one.gene$yday <- date.info$yday
  one.gene$season <- date.info$season
  one.gene$invasion <- date.info$invasion

  gene.info.list <- c(gene.info.list, list(one.gene))
}

gene.info.table <- rbindlist(l = gene.info.list)
rm(gene.info.list)

gene.info.table[ , genome := .(sub(pattern = "_scaffold_.*$",replacement = "", x = scaffold))]

# ---- save output ----

all.genomes <- unique(gene.info.table$genome)

for (g in all.genomes){
  one.genome <- gene.info.table[genome == g]
  fwrite(x = one.genome, file = paste0(g,"_gene_info_combined.tsv.gz"), sep = "\t", nThread = num.threads, compress = "gzip")
}



