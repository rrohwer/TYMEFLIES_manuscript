# RRR

# Example cov files from mapping year 2000 to ME2000-03-15pf_3300044539 on helheim
folder.path <- "data/2022-06-02_cov_files_for_ME2000-03-15pf_3300044539/covfiles/"
save.figures.here <- "figures/2022-06-02_coverage_clusters"
save.figures.here <- "~/Desktop/cov_clusts_1000"
my.files <- list.files(folder.path)
file.1 <- read.table(file = file.path(folder.path,my.files[1]), sep = "\t", comment.char = "", header = TRUE)
my.cov <- matrix(data = NA,nrow = nrow(file.1), ncol = length(my.files))
row.names(my.cov) <- file.1$X.ID  
colnames(my.cov) <- sub(pattern = "\\.cov", replacement = "", x = my.files)
my.cov[ ,1] <- file.1$Avg_fold
for (f in 2:length(my.files)){
  curr.file <- read.table(file = file.path(folder.path,my.files[f]), sep = "\t", comment.char = "", header = TRUE)
  cat(all.equal(row.names(my.cov), curr.file$X.ID),"\n")
  my.cov[ ,f] <- curr.file$Avg_fold
}
rm(curr.file, file.1)


# normalize by contig max abundance

contig.maxes <- apply(X = my.cov, MARGIN = 1, FUN = max)
index <- which(contig.maxes == 0) # 445 scaffolds have zero coverage- how??
contig.maxes <- contig.maxes[-index]
norm.cov <- my.cov[-index, ]
norm.cov <- norm.cov / contig.maxes * 100
head(norm.cov)


# cluster contigs using simplest kmeans algorithm

y <- kmeans(x = norm.cov, centers = 1000, iter.max = 1000, nstart = 1, algorithm = "Lloyd")

y$cluster[1:5] # which cluster each contig belongs to
head(y$centers) # matrix of cluster means- the average value for each cluster on each day
y$size # number of datapoints in each cluster
hist(y$size)

# look at clusters

clust = 24
for (clust in 1:length(y$size)){
  index <- which(y$cluster == clust)
  num.contigs <- y$size[clust]
  x.dates <- substr(x = colnames(my.cov), start = 3, stop = 12)
  x.dates <- parse_date_time(x = x.dates, orders = "ymd", tz = "Etc/GMT-5")
  png(filename = paste0(save.figures.here,"/cluster_",clust), width = 4, height = 4, units = "in", res = 100)
  par(mar = c(3.5,5,3,.5), oma = c(.1,.1,.1,.1))
  plot(x = c(min(x.dates),max(x.dates)), y = c(0,100), type = "n", ann = F, axes = F)
  for (r in 1:length(index)){
    lines(x = x.dates, y = norm.cov[index[r], ], col = adjustcolor("black",.1))
  }
  axis(side = 1, at = x.dates, labels = paste(month(x.dates),day(x.dates), sep = "-"), las = 2)
  axis(side = 2, las = 2)
  # mtext(text = "Sample", side = 1, line = 2.25)
  mtext(text = "Contig Coverage\n(% of contig max)", side = 2,line = 2.5)
  mtext(text = paste("Cluster",clust), side = 3, line = 1.5, cex = 1.2)
  mtext(text = paste(num.contigs,"Total Contigs"), side = 3, line = .5)
  mtext(text = "Year 2000,\nassembly 3/15/2000", side = 3, line = -2.5, outer = T, at = .01, adj = 0, cex = .7)
  dev.off()
}
