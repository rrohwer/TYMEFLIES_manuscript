# RRR
# approaching my time-series binning idea, try using the cov file data

library("vegan")

cov.a <- read.table(file = "data/2022-02-28_targeted_binning_testing/ME2015-06-10_3300042356.cov", 
                    sep = "\t", comment.char = "", header = TRUE)
cov.b <- read.table(file = "data/2022-02-28_targeted_binning_testing/ME2017-06-16_3300042330_to_ME2015-06-10_3300042356.cov",
                    sep = "\t", comment.char = "", header = TRUE)
cov.c <- read.table(file = "data/2022-02-28_targeted_binning_testing/ME2018-06-07_3300034013_to_ME2015-06-10_3300042356.cov",
                    sep = "\t", comment.char = "", header = TRUE)

# cov.a <- read.table(file = "ME2015-06-10_3300042356.cov", 
#                     sep = "\t", comment.char = "", header = TRUE)
# cov.b <- read.table(file = "ME2017-06-16_3300042330_to_ME2015-06-10_3300042356.cov",
#                     sep = "\t", comment.char = "", header = TRUE)
# cov.c <- read.table(file = "ME2018-06-07_3300034013_to_ME2015-06-10_3300042356.cov",
#                     sep = "\t", comment.char = "", header = TRUE)

colnames(cov.a)
all.equal(cov.a$X.ID, cov.b$X.ID)
all.equal(cov.a$X.ID, cov.c$X.ID)

my.cov <- cbind(cov.a$X.ID, cov.a$Avg_fold, cov.b$Avg_fold, cov.c$Avg_fold)
colnames(my.cov) <- c("contig","A","B","C")
temp <- my.cov[ ,1]
my.cov <- my.cov[ ,-1]
my.cov <- apply(X = my.cov, MARGIN = 2, FUN = as.numeric)
row.names(my.cov) <- temp
head(my.cov)

rm(cov.a, cov.b, cov.c, temp)

hist(my.cov)

# save example for Val
write.table(x = my.cov, file = "~/Desktop/example_coverage.tsv", quote = F, sep = "\t", row.names = T, col.names = T)

# didn't work
abund.cutoff <- 90
samples.cutoff <- 3
pres.abs <- my.cov >= abund.cutoff
head(pres.abs)
index.keep <- rowSums(pres.abs) >= samples.cutoff
trim.cov <- my.cov[index.keep, ]
dim(trim.cov)

my.dist <- dist(x = trim.cov, method = "euclidean", diag = FALSE, upper = FALSE)

my.dist <- vegdist(x = my.cov, method = "euclidean")

small.cov <- my.cov[1:10000, ]

my.dist <- dist(x = small.cov)

# this is not going to work because can't make the distance matrix, even on the server
# look into rRNA seq gene abundance analysis tools- they have a lot of genes too and are just basing it on abundance.

# did work

# https://stackoverflow.com/questions/21382681/kmeans-quick-transfer-stage-steps-exceeded-maximum
?kmeans()
x <- kmeans(x = my.cov, centers = 100, iter.max = 1000, nstart = 1, algorithm = "Lloyd")
head(my.cov)
dim(my.cov)
x$cluster[1:5] # which cluster each point belongs to
head(x$centers) # matrix of cluster means- the average value for each cluster on each day
dim(x$centers)
x$totss # total sum of squares- total variance
x$withinss # within-cluster sum of squares- within-cluster variance
x$tot.withinss # total sum of squares of all clusters
x$betweenss # total sum of squares between clusters
x$size # number of datapoints in each cluster
hist(x$size)

# normalize first
scaff.maxes <- apply(X = my.cov, MARGIN = 1, FUN = max)
summary(scaff.maxes)
index <- which(scaff.maxes == 0) # 445 scaffolds have zero coverage- how??
scaff.maxes <- scaff.maxes[-index]
norm.cov <- my.cov[-index, ]
norm.cov <- norm.cov / scaff.maxes * 100
head(norm.cov)
y <- kmeans(x = norm.cov, centers = 100, iter.max = 1000, nstart = 1, algorithm = "Lloyd")
head(norm.cov)
dim(norm.cov)
y$cluster[1:5] # which cluster each point belongs to
head(y$centers) # matrix of cluster means- the average value for each cluster on each day
dim(y$centers)
y$totss # total sum of squares- total variance
y$withinss # within-cluster sum of squares- within-cluster variance
y$tot.withinss # total sum of squares of all clusters
y$betweenss # total sum of squares between clusters
y$size # number of datapoints in each cluster
hist(y$size)

clust = 85
for (clust in 1:100){
  index <- which(y$cluster == clust)
  plot(x = c(1,3), y = c(0,100), type = "n")
  for (r in 1:length(index)){
    lines(x = 1:3, y = norm.cov[index[r], ], col = adjustcolor("black",.1))
  }
}

