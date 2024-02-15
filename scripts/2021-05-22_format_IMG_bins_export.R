# RRR
# I downloaded the bin stats from IMG, so take big picture looks at them
# I know their quality, the sample they're from, and the GTDB taxonomy

# ---- set-up and import ----

library(readxl)
library(tidyr)

clean.rds.output <- "data/2021-05-22_bin_stats/all_img_bins.rds"
clean.tsv.output <- "data/2021-05-22_bin_stats/all_img_bins.tsv"

ids.only.output.allbins <- "data/2021-05-22_bin_stats/IDs_only_all_bins.txt"
ids.only.output.HQ <- "data/2021-05-22_bin_stats/IDs_only_HQ_bins.txt"
ids.only.output.MQ <- "data/2021-05-22_bin_stats/IDs_only_MQ_bins.txt"

# downloaded separately b/c couldn't remove the 5000 row limit
IMG_HQ <- read_excel(path = "data/2021-05-22_bin_stats/2021-05-20_onlyHQ_BinIDs_for_504350.xlsx")
IMG_MQ1 <- read_excel(path = "data/2021-05-22_bin_stats/2021-05-22_MQbins_scaff_0-200.xlsx")
IMG_MQ2 <- read_excel(path = "data/2021-05-22_bin_stats/2021-05-22_MQbins_scaff_200-300.xlsx")
IMG_MQ3 <- read_excel(path = "data/2021-05-22_bin_stats/2021-05-22_MQbins_scaff_300-3280.xlsx")

IMG_all <- rbind(IMG_HQ, IMG_MQ1, IMG_MQ2, IMG_MQ3)

# ---- parse out my meaningful sample names ----

my.IDs <- IMG_all[ ,3, drop =T] %>%
  sub(pattern = "Freshwater microbial communities from Lake Mendota, Madison, Wisconsin, United States - TYMEFLIES-", replacement = "", x = .) %>%
  sub(pattern = "Freshwater microbial communities from Lake Mendota, Wisconsin, United States - TYMEFLIES-", replacement = "", x = .)

Sample.Name <- sub(pattern = "-.*$", replacement = "", x = my.IDs)

TF.ex.Code <- sub(pattern = "^.*-", replacement = "", x = my.IDs)
rr <- substr(x = TF.ex.Code, start = 1, stop = 2)
ex.num <- substr(x = TF.ex.Code, start = 3, stop = 7) %>%
  as.numeric(.) # these NAs are gen donors, take care of in next code block
Extraction.Code <- paste(rr, ex.num)

index <- which(Sample.Name == "CONTROL")
Sample.Name[index] <- "ME08Nov2018GD"
Extraction.Code[index] <- TF.ex.Code[index] %>%
  sub(pattern = "^GENDONOR", replacement = "gd", x = .)

cbind(1:ncol(IMG_all), colnames(IMG_all))
IMG_all[1:3,20]

bins <- data.frame("Bin.ID" = IMG_all[ ,2,drop=T],
                   "IMG.Taxon.ID" = as.character(IMG_all[ ,4,drop=T]),
                   "ITS.Proposal.ID" = as.character(IMG_all[ ,20,drop=T]),
                   "Sample.Name" = Sample.Name,
                   "Extraction.Code" = Extraction.Code,
                   "Bin.Quality" = IMG_all[ ,5,drop=T],
                   "Completeness" = IMG_all[ ,11,drop=T],
                   "Contamination" = IMG_all[ ,12,drop=T],
                   "Number.Bases" = IMG_all[ ,13,drop=T],
                   "Number.Scaffolds" = IMG_all[ ,19,drop=T],
                   "Number.Genes" = IMG_all[ ,18,drop=T],
                   "rRNA.5S" = IMG_all[ ,14,drop=T],
                   "rRNA.16S" = IMG_all[ ,15,drop=T],
                   "rRNA.23S" = IMG_all[ ,16,drop=T],
                   "tRNA.genes" = IMG_all[ ,17,drop=T])

# ---- fix up taxonomy names ----

# stealing this function from the taxass script reformat_taxonomy_nomenclature.R
make.degenerates.unique <- function(tax, voldemort){
  # tax is the taxonomy table without seqid column
  # voldemort is whatever text you're making unique (b/c it cannot be named... )
  for (t in 1:ncol(tax)){
    index <- which(tax[ ,t] == voldemort)
    tax[index,t] <- paste(tax[index,t], tax[index,t - 1], sep = ".")
  }
  
  remove.extra.u <- function(x){
    dot.voldemort <- paste0(".", voldemort)
    x <- gsub(pattern = dot.voldemort, replacement = "", x = x)
    return(x)
  }
  tax <- apply(X = tax, MARGIN = 2, FUN = remove.extra.u)
  
  return(tax)
}

# oh the other reformat functions are only for unique taxonomy tables- this has repeat orgs so repeat rows.

gtdb <- IMG_all[ ,7]
gtdb <- separate(data = gtdb, col = 1, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = "; ", remove = TRUE, fill = "right", extra = "warn")
gtdb <- as.matrix(gtdb)
index <- is.na(gtdb)
gtdb[index] <- "unclassified"
gtdb <- make.degenerates.unique(tax = gtdb, voldemort = "unclassified")

img.tax <- IMG_all[ ,6]
img.tax <- separate(data = img.tax, col = 1, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = "; ", remove = TRUE, fill = "right", extra = "warn")
img.tax <- as.matrix(img.tax)
index <- is.na(img.tax)
img.tax[index] <- "unclassified"
img.tax <- make.degenerates.unique(tax = img.tax, voldemort = "unclassified")

colnames(gtdb) <- paste0("GTDB.",colnames(gtdb))
colnames(img.tax) <- paste0("IMG.", colnames(img.tax))

bins <- cbind(bins, gtdb, img.tax)

str(bins)

all.IDs <- bins$Bin.ID
index <- bins$Bin.Quality == "HQ"
HQ.IDs <- bins$Bin.ID[index]
MQ.IDs <- bins$Bin.ID[!index]
unique(bins$Bin.Quality[index])
unique(bins$Bin.Quality[!index])

# ---- export formatted bin stats ----

# saveRDS(object = bins, file = clean.rds.output)
# write.table(x = bins, file = clean.tsv.output, quote = FALSE, sep = "\t", row.names = FALSE)
# 
# write.table(x = all.IDs, file = ids.only.output.allbins, quote = F, row.names = F)
# write.table(x = HQ.IDs, file = ids.only.output.HQ, quote = F, row.names = F)
# write.table(x = MQ.IDs, file = ids.only.output.MQ, quote = F, row.names = F)




# ~end~