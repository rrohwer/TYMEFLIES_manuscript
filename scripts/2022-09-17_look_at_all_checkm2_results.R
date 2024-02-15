# RRR

data.folder <- "data/2022-09-17_checkm2_output/checkm2_quality_files/"

combo <- NULL
for(f in list.files(data.folder)){
  my.tab <- read.delim(file = file.path(data.folder,f))
  combo <- rbind(combo, my.tab)
}

hist(combo$Completeness, breaks = 100)
hist(combo$Contamination, breaks = 1000)
# zoom in
hist(combo$Contamination[combo$Contamination <= 20], breaks = 40)

summary(combo$Completeness)
summary(combo$Contamination)


num.bins <- matrix(data = NA, nrow = 4, ncol = 10)
completeness <- seq.int(from = 10, to = 100, by = 10)
redundancy <- c(1,2.5,5,10)
colnames(num.bins) <- paste("over",completeness, sep = "_")
row.names(num.bins) <- paste("under", redundancy, sep = "_")

for (c in 1:length(completeness)){
  for (r in 1:length(redundancy)){
    num.bins[r,c] <- sum(combo$Completeness >= completeness[c] & combo$Contamination <= redundancy[r])
  }
}


tot.HQ <- sum(combo$Completeness >= 50 & combo$Contamination <= 10)

jgi <- readRDS("data/2021-05-22_bin_stats/all_img_bins.rds") # only 12,231

num.bins

saveRDS(object = combo, file = "data/2022-09-17_checkm2_output/combined_checkm2_quality_files.rds")

write.csv(x = num.bins, file = "data/2022-09-17_checkm2_output/num_bins.csv", row.names = T)
