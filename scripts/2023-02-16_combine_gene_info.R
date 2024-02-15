# RRR
# read in the gene info tables from instrain
# combine them into a single file and save
# edited the genome info script to go faster with data.table fread and rbindlist

library(lubridate)
library(limony)
data("seasons")
library(data.table)

folder.path <- "data/2023-02-12_instrain_results_TIMEOUT"

my.files <- list.files(path = file.path(folder.path,"gene_info_files"), full.names = T)

get.date.info <- function(sample.name, seasons, key){
  sample.name <- sub(pattern = "^.*gene_info_files/", replacement = "", x = sample.name)
  sample.name <- sub(pattern = ".IS_gene_info.tsv", replacement = "", x = sample.name)
  
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
  
  # sample.date <- as.character(sample.date) # for faster fwrite
  
  return(list("name" = sample.name,"date" = sample.date, 
              "season" = my.season, "color.season" = my.season.color, 
              "invasion" = invasion, "color.invasion" = invasion.color))
}

gene.info.list <- list()
for (f in my.files){
  one.gene <- fread(file = f, sep = "\t", header = T)
  date.info <- get.date.info(sample.name = f, seasons = seasons, key = key)
  one.gene$sample <- date.info$name
  one.gene$date <- date.info$date
  one.gene$season <- date.info$season
  one.gene$invasion <- date.info$invasion
  one.gene$color.season <- date.info$color.season
  one.gene$color.invasion <- date.info$color.invasion
  
  gene.info.list <- c(gene.info.list, list(one.gene))
}

gene.info.table <- rbindlist(l = gene.info.list)
rm(gene.info.list)

gene.info.table[ , genome := .(sub(pattern = "_scaffold_.*$",replacement = "", x = scaffold))]

# fwrite(x = gene.info.table, file = file.path(folder.path,"gene_info_combined.tsv.gz"), sep = "\t", compress = "gzip")
saveRDS(object = gene.info.table, file = file.path(folder.path,"gene_info_combined.rds"))
