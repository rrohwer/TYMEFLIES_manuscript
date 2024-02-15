# RRR
library(lubridate)
library(readxl)

bigtyme <- read_excel("data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx")

output.folder <- "data/2022-06-27_binning_groups/"

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

# ---- parse bigtyme sample info ----

bigtyme <- bigtyme[!is.na(bigtyme$taxon_oid), ]
bigtyme$sample.names <- paste0(bigtyme$TYMEFLIES.name, "_",bigtyme$taxon_oid)
bigtyme$sample.dates <- parse_date_time(x = paste(bigtyme$Year,bigtyme$Month,bigtyme$Day), orders = "ymd", tz = "Etc/GMT-5")
bigtyme <- bigtyme[order(bigtyme$sample.dates), ]

for (y in unique(bigtyme$Year)){
  cat(y, " ", sum(bigtyme$Year == y),"\n")
}

# ---- make a list of sample groups to cross-map ----

# aim for ~50 samples per group, go by year
# split these out manually

mapping.groups <- make.empty.list.structure(c("test",paste0("group",1:8)))
bigtyme$mapping.group <- NA

index <- which(bigtyme$Year == 2000)
mapping.groups$test <- bigtyme$sample.names[index]

index <- which(bigtyme$Year == 2000 | bigtyme$Year == 2001 | bigtyme$Year == 2002)
mapping.groups$group1 <- bigtyme$sample.names[index]
bigtyme$mapping.group[index] <- "group1"

index <- which(bigtyme$Year == 2003 | bigtyme$Year == 2004 | bigtyme$Year == 2005 | bigtyme$Year == 2006 | bigtyme$Year == 2007)
mapping.groups$group2 <- bigtyme$sample.names[index]
bigtyme$mapping.group[index] <- "group2"

index <- which(bigtyme$Year == 2008 | bigtyme$Year == 2009 | bigtyme$Year == 2010 | bigtyme$Year == 2011)
mapping.groups$group3 <- bigtyme$sample.names[index]
bigtyme$mapping.group[index] <- "group3"

index <- which(bigtyme$Year == 2012 | bigtyme$Year == 2013)
mapping.groups$group4 <- bigtyme$sample.names[index]
bigtyme$mapping.group[index] <- "group4"

index <- which(bigtyme$Year ==  2014)
mapping.groups$group5 <- bigtyme$sample.names[index]
bigtyme$mapping.group[index] <- "group5"

index <- which(bigtyme$Year == 2015)
mapping.groups$group6 <- bigtyme$sample.names[index]
bigtyme$mapping.group[index] <- "group6"

index <- which(bigtyme$Year == 2016 | bigtyme$Year == 2017)
mapping.groups$group7 <- bigtyme$sample.names[index]
bigtyme$mapping.group[index] <- "group7"

index <- which(bigtyme$Year == 2018 | bigtyme$Year == 2019)
mapping.groups$group8 <- bigtyme$sample.names[index]
bigtyme$mapping.group[index] <- "group8"

# ---- estimate timing ----

num.samples <- NULL
for (e in names(mapping.groups)){
  num.samples <- c(num.samples, length(mapping.groups[[e]]))
}
names(num.samples) <- names(mapping.groups)
num.samples

node.hours <- num.samples ^ 2 * .4
sum(node.hours)

# ---- make job submission checklist ----

colnames(bigtyme)
checklist <- cbind(bigtyme[ ,c("sample.names","mapping.group")], 
                   "mapping.submitted" = "", "completed" = "", "N" = "","n" = "","my.time" = "","node.hours" = "", "notes.mapping" ="",
                   "sorting.submitted" = "", "completed" = "", "N" = "","n" = "","my.time" = "","node.hours" = "", "notes.sorting" ="", 
                   "deleted.sams" = "",
                   "depths.submitted" = "", "completed" = "", "N" = "","n" = "","my.time" = "","node.hours" = "", "notes.depths" ="",
                   "metabat.submitted" = "", "completed" = "", "N" = "","n" = "","my.time" = "","node.hours" = "", "notes.metabat" ="",
                   "concoct.submitted" = "", "completed" = "", "N" = "","n" = "","my.time" = "","node.hours" = "", "notes.concoct" ="",
                   "maxbin.submitted" = "", "completed" = "", "N" = "","n" = "","my.time" = "","node.hours" = "", "notes.maxbin" ="",
                   "dastool.submitted" = "", "completed" = "", "N" = "","n" = "","my.time" = "","node.hours" = "", "notes.depths" ="",
                   "deleted.bams" = "")

# ---- export ----

saveRDS(object = bigtyme, file = paste0(output.folder,"2022-06-27_bigtyme.rds"))

saveRDS(object = mapping.groups, file = paste0(output.folder, "mapping_groups_list.rds"))

write.csv(x = checklist, file = paste0(output.folder, "1_mapping_checklist.csv"), row.names = F)
