# RRR

save.figures.here <- "figures/2022-05-11_coverage_clusters/"

# import coverage files

cov.a <- read.table(file = "data/2022-02-28_targeted_binning_testing/ME2015-06-10_3300042356.cov", sep = "\t", comment.char = "", header = TRUE)
cov.b <- read.table(file = "data/2022-02-28_targeted_binning_testing/ME2017-06-16_3300042330_to_ME2015-06-10_3300042356.cov", sep = "\t", comment.char = "", header = TRUE)
cov.c <- read.table(file = "data/2022-02-28_targeted_binning_testing/ME2018-06-07_3300034013_to_ME2015-06-10_3300042356.cov", sep = "\t", comment.char = "", header = TRUE)


# merge coverage files

colnames(cov.a)
all.equal(cov.a$X.ID, cov.b$X.ID)
all.equal(cov.a$X.ID, cov.c$X.ID)

my.cov <- cbind(cov.a$X.ID, cov.a$Avg_fold, cov.b$Avg_fold, cov.c$Avg_fold)
my.cov <- matrix(c(cov.a$Avg_fold, cov.b$Avg_fold, cov.c$Avg_fold), ncol = 3)
rownames(my.cov) <- cov.a$X.ID
colnames(my.cov) <- c("A","B","C")
head(my.cov)

rm(cov.a, cov.b, cov.c)


# normalize by contig max abundance

contig.maxes <- apply(X = my.cov, MARGIN = 1, FUN = max)
index <- which(contig.maxes == 0) # 445 scaffolds have zero coverage- how??
contig.maxes <- contig.maxes[-index]
norm.cov <- my.cov[-index, ]
norm.cov <- norm.cov / contig.maxes * 100
head(norm.cov)


# cluster using simplest kmeans algorithm

y <- kmeans(x = norm.cov, centers = 100, iter.max = 1000, nstart = 1, algorithm = "Lloyd")

y$cluster[1:5] # which cluster each point belongs to
head(y$centers) # matrix of cluster means- the average value for each cluster on each day
y$size # number of datapoints in each cluster


# look at clusters

clust = 24
for (clust in 1:100){
  index <- which(y$cluster == clust)
  num.contigs <- y$size[clust]
  png(filename = paste0(save.figures.here,"/cluster_",clust), width = 4, height = 4, units = "in", res = 100)
  par(mar = c(3.5,5.5,4,1), oma = c(.1,.1,.1,.1))
  plot(x = c(1,3), y = c(0,100), type = "n", ann = F, axes = F)
  for (r in 1:length(index)){
    lines(x = 1:3, y = norm.cov[index[r], ], col = adjustcolor("black",.1))
  }
  axis(side = 1, at = c(1:3), labels = c("A","B","C"))
  axis(side = 2, las = 2)
  mtext(text = "Sample", side = 1, line = 2.25)
  mtext(text = "Contig Coverage\n(% of contig max)", side = 2,line = 3)
  mtext(text = paste("Cluster",clust), side = 3, line = 2, cex = 1.2)
  mtext(text = paste(num.contigs,"Total Contigs"), side = 3, line = 1)
  dev.off()
}
