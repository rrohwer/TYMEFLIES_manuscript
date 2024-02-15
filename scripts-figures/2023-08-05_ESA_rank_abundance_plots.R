# RRR
# just how abundant are my "cherry-picked" genomes?

library(data.table)
library(lubridate)

coverm <- fread(input = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")


colnames(coverm)[3] <- "rel.abund"
coverm <- coverm[ ,1:3]
coverm[ ,date := substr(Sample, start = 3, stop = 12)]
coverm[ ,date := parse_date_time(x = date, orders = "ymd")]
coverm <- coverm[Genome != "unmapped"]

coverm <- dcast(data = coverm, formula = Genome ~ date, fun.aggregate = mean, value.var = "rel.abund")
rank.abund <- data.table("Genome" = coverm$Genome, 
                     "Sum.Abund" = rowSums(coverm[ ,-1]), 
                     "Mean.Abund" = apply(X = coverm[ ,-1], MARGIN = 1, FUN = mean), 
                     "Med.Abund" = apply(X = coverm[ ,-1], MARGIN = 1, FUN = median))
rank.abund <- rank.abund[order(Mean.Abund, decreasing = T)]

# just add abundances to drep table, will be good in the supplement

