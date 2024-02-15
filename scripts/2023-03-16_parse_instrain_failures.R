# RRR
# grep out file sizes and job numbers and job times

jobnames <- read.delim(file = "data/2023-03-13_inStrain_on_drep96/errors/job_names.txt", header = F, sep = "\t")
jobtimes <- read.delim(file = "data/2023-03-13_inStrain_on_drep96/errors/job_completion_times.txt", header = F, sep = "\t")
bamsizes <- read.delim(file = "data/2023-03-13_inStrain_on_drep96/errors/bam_sizes.txt", header = F, sep = "\t")

head(jobnames)
head(jobtimes)
x <- merge(x = jobnames, y = jobtimes, by.x = "V2", by.y = "V2", all = T)
head(x)
x <- x[ ,c(3,1,5)]
colnames(x) <- c("sample","job","seconds")
colnames(bamsizes) <- c("GB","sample")
head(bamsizes)
x <- merge(x = x, y = bamsizes, by = "sample", all = T)
head(x)
x$hours <- x$seconds/60/60
x$GB <- sub(pattern = "G", replacement = "", x = x$GB) |>
  as.numeric()

plot(x = x$GB, y = x$hours, xlab = "Input file size (GB)", ylab = "Job completion time (Hours)", main = "InStrain 60-node TIMEOUT run")
points(x = x$GB[x$sample == "ME2018-06-14_3300042564"], y = x$hours[x$sample == "ME2018-06-14_3300042564"], pch = 19, col = "red")
text(x = x$GB[x$sample == "ME2018-06-14_3300042564"], y = x$hours[x$sample == "ME2018-06-14_3300042564"] + .4, labels = round(x$hours[x$sample == "ME2018-06-14_3300042564"],1), col = "red")
# took 3 hrs 30 min in test run on 96 threads

x$sample[which(x$hours == (max(x$hours, na.rm = T)))]
x$GB[which(x$hours == (max(x$hours, na.rm = T)))]
# slowest one to complete was "ME2017-07-24_3300042509", 9.7 GB bam file, 9.5 hrs

plot(x = x$GB, y = x$hours, xlab = "Input file size (GB)", ylab = "Job completion time (Hours)", main = "InStrain 60-node TIMEOUT run", type = "n")
text(x = x$GB, y = x$hours, labels = x$job, cex = .5)
plot(x = x$GB, y = x$job, xlab = "Input file size (GB)", ylab = "Job number", main = "InStrain 60-node TIMEOUT run")

