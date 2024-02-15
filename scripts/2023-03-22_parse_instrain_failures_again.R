# RRR
# grep out file sizes and job numbers and job times

# ---- does bam file size determine timing? ----

bam.sizes <- read.delim(file = "data/2023-03-13_inStrain_on_drep96/errors/bam_sizes.txt", header = F, sep = "\t")
colnames(bam.sizes) <- c("GB","Sample")
bam.sizes$GB <- sub("G","",bam.sizes$GB) |>
  as.numeric()

# ---- main instrain submission (60 nodes) ----

completed <- read.table("data/2023-03-13_inStrain_on_drep96/errors/Job_completed_list.txt", sep = " ")
commands <- read.table("data/2023-03-13_inStrain_on_drep96/errors/job_commands_list.txt", sep = " ", skip = 1)

head(completed)
completed <- completed[ ,c(3,6,4)]
colnames(completed) <- c("Job","seconds","status")
completed$hours <- completed$seconds / 60 / 60
unique(completed$status) # these 395 all finished

head(commands)
commands <- commands[ ,c(3,6,21)]
colnames(commands) <- c("Task","Job","Sample")
commands$Sample <- sub("\\.IS", "",commands$Sample)

completed <- merge(x = completed, y = commands, by = "Job", all.x = TRUE, all.y = FALSE)
completed <- merge(x = completed, y = bam.sizes, by = "Sample", all.x = TRUE, all.y = FALSE)
failure <- completed
head(failure)

# ---- repeat unfinished instrain submission (60 nodes) ----

completed <- read.table("data/2023-03-13_inStrain_on_drep96/errors_unfinished/Job_complete_list_unfinished.txt", sep = " ")
commands <- read.table("data/2023-03-13_inStrain_on_drep96/errors_unfinished/job_commands_list_unfinished.txt", sep = " ", skip = 1)

head(completed)
completed <- completed[ ,c(3,6,4)]
colnames(completed) <- c("Job","seconds","status")
completed$hours <- completed$seconds / 60 / 60
unique(completed$status) # these 57 all finished

head(commands)
commands <- commands[ ,c(3,6,21)]
colnames(commands) <- c("Task","Job","Sample")
commands$Sample <- sub("\\.IS", "",commands$Sample)
nrow(commands) # these 67 were launched

completed <- merge(x = completed, y = commands, by = "Job", all.x = TRUE, all.y = FALSE)
completed <- merge(x = completed, y = bam.sizes, by = "Sample", all.x = TRUE, all.y = FALSE)
failure.1 <- completed
head(failure.1)

# ---- look at how long runs that finished took ----

png(filename = "figures/2023-03-22_instrain_job_completion_times/Time_vs_BamSize.png", width = 6.5, height = 4, units = "in", res = 300)

par(mfrow = c(1,2), mar = c(4,4,2,.5), oma = c(.1,.1,3,.1))

y.lim <- c(1,9.75)

plot(x = failure$GB, y = failure$hours, xlab = "Input file size (GB)", ylab = "Job completion time (Hours)", ylim = y.lim)
points(x = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y = failure$hours[failure$Sample == "ME2018-06-14_3300042564"], pch = 19, col = "red")
# took 3:30:39 min in test run 
largest.hrs <- 3 + 30 / 60 + 39 / 60 / 60
points(x = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y = largest.hrs, pch = 19, col = "blue")
segments(x0 = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], x1 = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y0 = largest.hrs, y1 = failure$hours[failure$Sample == "ME2018-06-14_3300042564"], col = "red")
points(x = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y = failure$hours[failure$Sample == "ME2017-07-24_3300042509"], pch = 19, col = "red")
# took 03:39:45 min in test run 
slowest.hrs <- 3 + 39 / 60 + 45 / 60 / 60
points(x = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y = largest.hrs, pch = 19, col = "blue")
segments(x0 = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], x1 = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y0 = largest.hrs, y1 = failure$hours[failure$Sample == "ME2017-07-24_3300042509"], col = "red")
mtext(text = "InStrain 60-node TIMEOUT runs", outer = T, line = 1)
mtext(text = "Run all", outer = F, line = 1)

plot(x = failure.1$GB, y = failure.1$hours, xlab = "Input file size (GB)", ylab = "Job completion time (Hours)", ylim = y.lim)
mtext(text = "Repeat unfinished", outer = F, line = 1)

dev.off()

# ----

png(filename = "figures/2023-03-22_instrain_job_completion_times/Time_vs_JobNumber.png", width = 6.5, height = 4, units = "in", res = 300)

par(mfrow = c(1,2), mar = c(4,4,2,.5), oma = c(.1,.1,3,.1))

y.lim <- c(1,9.75)

plot(x = failure$Job, y = failure$hours, xlab = "Job Number", ylab = "Job completion time (Hours)", ylim = y.lim)
points(x = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y = failure$hours[failure$Sample == "ME2018-06-14_3300042564"], pch = 19, col = "red")
# took 3:30:39 min in test run 
largest.hrs <- 3 + 30 / 60 + 39 / 60 / 60
points(x = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y = largest.hrs, pch = 19, col = "blue")
segments(x0 = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], x1 = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y0 = largest.hrs, y1 = failure$hours[failure$Sample == "ME2018-06-14_3300042564"], col = "red")
points(x = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y = failure$hours[failure$Sample == "ME2017-07-24_3300042509"], pch = 19, col = "red")
# took 03:39:45 min in test run 
slowest.hrs <- 3 + 39 / 60 + 45 / 60 / 60
points(x = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y = largest.hrs, pch = 19, col = "blue")
segments(x0 = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], x1 = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y0 = largest.hrs, y1 = failure$hours[failure$Sample == "ME2017-07-24_3300042509"], col = "red")
mtext(text = "InStrain 60-node TIMEOUT runs", outer = T, line = 1)
mtext(text = "Run all", outer = F, line = 1)

plot(x = failure.1$Job, y = failure.1$hours, xlab = "Job Number", ylab = "Job completion time (Hours)", ylim = y.lim)
mtext(text = "Repeat unfinished", outer = F, line = 1)

dev.off()

# ----

png(filename = "figures/2023-03-22_instrain_job_completion_times/Time_vs_TaskNumber.png", width = 6.5, height = 4, units = "in", res = 300)

par(mfrow = c(1,2), mar = c(4,4,2,.5), oma = c(.1,.1,3,.1))

y.lim <- c(1,9.75)

plot(x = failure$Task, y = failure$hours, xlab = "Task Number", ylab = "Job completion time (Hours)", ylim = y.lim)
points(x = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y = failure$hours[failure$Sample == "ME2018-06-14_3300042564"], pch = 19, col = "red")
# took 3:30:39 min in test run 
largest.hrs <- 3 + 30 / 60 + 39 / 60 / 60
points(x = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y = largest.hrs, pch = 19, col = "blue")
segments(x0 = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], x1 = failure$GB[failure$Sample == "ME2018-06-14_3300042564"], y0 = largest.hrs, y1 = failure$hours[failure$Sample == "ME2018-06-14_3300042564"], col = "red")
points(x = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y = failure$hours[failure$Sample == "ME2017-07-24_3300042509"], pch = 19, col = "red")
# took 03:39:45 min in test run 
slowest.hrs <- 3 + 39 / 60 + 45 / 60 / 60
points(x = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y = largest.hrs, pch = 19, col = "blue")
segments(x0 = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], x1 = failure$GB[failure$Sample == "ME2017-07-24_3300042509"], y0 = largest.hrs, y1 = failure$hours[failure$Sample == "ME2017-07-24_3300042509"], col = "red")
mtext(text = "InStrain 60-node TIMEOUT runs", outer = T, line = 1)
mtext(text = "Run all", outer = F, line = 1)

plot(x = failure.1$Task, y = failure.1$hours, xlab = "Task Number", ylab = "Job completion time (Hours)", ylim = y.lim)
mtext(text = "Repeat unfinished", outer = F, line = 1)

dev.off()
