# RRR
# process the SNVs files into per-genome files that include all dates for that genome
# unlike gene_info, these are even bigger so don't load all into memory at once
# instead, load one at a time, write per-genome information to files, but append to the same files as go.
# need to run linearly since appending, but can fread with multiple threads to speed up.
# example:
# Rscript process_instrain/combine_SNP_info.R ../yggshare/current_projects/TYMEFLIES/tymeflies/runinstrain96 process_instrain/drep_results_all_genomes_0.96.rds process_instrain/limony_seasons.rds 30 per-genome_SNVs

library(lubridate)
library(data.table)

user.input <- commandArgs(trailingOnly = TRUE)
folder.path <- user.input[1]    # ../yggshare/current_projects/TYMEFLIES/tymeflies/runinstrain96/
sample.key <- user.input[2]     # drep_results_all_genomes_0.96.rds
season.key <- user.input[3]     # limony_seasons.rds
num.threads <- as.numeric(user.input[4])    # 40
output.folder <- user.input[5] # ./per-genome_SNVs

# # local test
# cat("\n\nyou forgot to comment out local paths!\n\n")
# folder.path <- "data/2023-07-16_SNVs_example_files/from_inStrain"
# sample.key <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"
# # library(limony)
# # data("seasons")
# # saveRDS(object = seasons, file = "data/2023-07-16_SNVs_example_files/limony_seasons.rds")
# season.key <- "data/2023-07-16_SNVs_example_files/limony_seasons.rds"
# num.threads <- 1
# output.folder <- "data/2023-07-16_SNVs_example_files/generated_per-genome"


# ---- functions ----

get.date.info <- function(sample.name, seasons){
  # get sample date and associated info for the gene_info file
  
  sample.name <- sub(pattern = "^.*output/*", replacement = "", x = sample.name)
  sample.name <- sub(pattern = ".IS_SNVs.*$", replacement = "", x = sample.name)
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

# ---- parse files ----
if (length(list.files(path = output.folder)) > 0){
  cat("\n\nSTOP! This script appends files, so the output folder must be EMPTY and it is not!\n\n")
}else{
  
  seasons <- readRDS(season.key)
  
  sample.names <- readRDS(sample.key)
  sample.names <- unique(sample.names$tymeflies.name)
  
  my.files <- list.files(path = paste0(folder.path,"/",sample.names,".IS/","output/"), 
                         pattern = "SNVs", full.names = TRUE) 
  
  for (f in my.files){
    # read in one SNVs file at a time
    cat("processing ",f,"\n")
    one.SNVs <- fread(file = f, sep = "\t", header = T, nThread = num.threads)
    
    # add date info for which metagenome was mapped
    date.info <- get.date.info(sample.name = f, seasons = seasons)
    one.SNVs$sample <- date.info$name
    one.SNVs$date <- date.info$date
    one.SNVs$year <- date.info$year
    one.SNVs$yday <- date.info$yday
    one.SNVs$season <- date.info$season
    one.SNVs$invasion <- date.info$invasion
    
    # add genome info for which genome was mapped to
    one.SNVs[ , genome := .(sub(pattern = "_scaffold_.*$",replacement = "", x = scaffold))]
    
    # export files, appending to each per-genome file as go
    for (g in unique(one.SNVs$genome)){
      per.genome <- one.SNVs[genome == g]
      fwrite(x = per.genome, file = file.path(output.folder,paste0(g,"_SNVs.tsv.gz")), sep = "\t", nThread = num.threads, compress = "gzip", append = TRUE)
    }
  }
  
}





