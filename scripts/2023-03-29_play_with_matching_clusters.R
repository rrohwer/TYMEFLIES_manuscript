# RRR
# goal is to classify as matches or lags
# correlate div vs. abund cluster patterns
# plot to manually look at correlation values and determine a cutoff for match or lag
# apply cutoff to classify clusters
# expand classification to each genome

library(data.table)

div.genomes <- fread(file = "data/2023-03-21_clustering_analysis_results/div_seasonal_table.tsv.gz")
abund.genomes <- fread(file = "data/2023-03-21_clustering_analysis_results/abund_seasonal_table.tsv.gz")

div.clusts <- div.genomes[ ,.(nuc.div = mean(nuc.div.perc)), by = .(div.cluster = cluster, season)]
abund.clusts <- abund.genomes[ ,.(abund.perc = mean(abund.perc.max)), by = .(abund.cluster = cluster, season)]

div.clusts[ , div.cluster := paste0("div.",div.cluster)]
abund.clusts[ , abund.cluster := paste0("abund.",abund.cluster)]

div.clusts <- dcast(data = div.clusts, formula = div.cluster ~ season, value.var = "nuc.div")
abund.clusts <- dcast(data = abund.clusts, formula = abund.cluster ~ season, value.var = "abund.perc")

corr.clusts <- data.table(clust.d = rep(div.clusts$div.cluster, each = length(abund.clusts$abund.cluster)),
                          clust.a = rep(abund.clusts$abund.cluster, times = length(div.clusts$div.cluster)),
                          cor = 0)

div.clusts <- as.matrix(x = div.clusts, rownames = T)
abund.clusts <- as.matrix(x = abund.clusts, rownames = T)

# pearson corr wants normal data
# if norm by max, all vectors have value 100 and it looks like a normal curve cut in the middle
# normalize cluster averages by mean cluster average instead of max to get a more normal distribution
# note (clustering done with max-based normalization, but Ward.D method does not expect normality)
div.clusts <- div.clusts / apply(X = div.clusts, MARGIN = 1, FUN = mean)
abund.clusts <- abund.clusts / apply(X = abund.clusts, MARGIN = 1, FUN = mean)

# for(r in 1:nrow(div.clusts)){
#   hist(x = div.clusts[r, ], breaks = 6)
# }
# for(r in 1:nrow(abund.clusts)){
#   hist(x = abund.clusts[r, ], breaks = 6)
# }
hist(div.clusts) 
hist(abund.clusts)

# can't figure out how to do without a loop, but not a ton of clusters so OK
for (r in 1:nrow(corr.clusts)){
  my.d <- corr.clusts[r, clust.d]
  my.a <- corr.clusts[r, clust.a]
  my.corr <- cor(x = div.clusts[row.names(div.clusts) == my.d, ], y = abund.clusts[row.names(abund.clusts) == my.a, ])
  corr.clusts[r,"cor"] <- my.corr
}

stripchart(x = corr.clusts$cor, method = "jitter", vertical = T, pch = 21)
hist(corr.clusts$cor)
# pretty evenly distributed btwn -1 to 1
# need to choose a cutoff for match vs lag
# plot to see at what level it looks clear:

for(r in 1:nrow(corr.clusts)){
  my.d <- div.clusts[row.names(div.clusts) == corr.clusts[r,clust.d], ]
  my.a <- abund.clusts[row.names(abund.clusts) == corr.clusts[r,clust.a], ]
  my.corr <- corr.clusts[r,cor]
  plot(x = c(1,6), y = c(min(c(my.a, my.d)), max(c(my.a,my.d))), type = "n", xlab = "season", ylab = "cluster values", main = round(my.corr,2))
  lines(x = 1:6, y = my.d, col = "blue", lwd = 3)
  lines(x = 1:6, y = my.a, col = "red3", lwd = 3)
}



