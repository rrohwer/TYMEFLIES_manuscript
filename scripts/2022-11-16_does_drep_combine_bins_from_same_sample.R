# RRR
# does drep ever combine bins from the same sample?
# would mean A. that it's probably too coarse
# and that B. concatenating all same-sample bins to map will still split bins

mybins <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")

mybins <- mybins[!is.na(mybins$drep.cluster.95ANI), ]

length(unique(mybins$drep.cluster.95ANI))

mytab <- data.frame("cluster" = unique(mybins$drep.cluster.95ANI), "num.bins" = NA, "num.samples" = NA)
for (clust in mytab$cluster){
  onebin <- mybins[mybins$drep.cluster.95ANI == clust, ]
  mytab[mytab$cluster == clust, "num.bins"] <- nrow(onebin)
  mytab[mytab$cluster == clust, "num.samples"] <- length(unique(onebin$tymeflies.name))
}
head(mytab)

index <- which(mytab$num.bins != mytab$num.samples)
index
mytab[index, ]

# one cluster contains two bins from the same sample

mybins[mybins$drep.cluster.95ANI == "135_1", ]

# so 95% was "too coarse" for one of the bacteroidota