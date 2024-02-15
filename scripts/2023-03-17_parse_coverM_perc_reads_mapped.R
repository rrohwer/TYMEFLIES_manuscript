# RRR

# do this to the output files:
# output $ cat coverM_drep96_terminal_output_minANI0.txt | grep sample > coverM_reads_mapped_0.txt
# output $ cat coverM_drep96_terminal_output_minANI93.txt | grep sample > coverM_reads_mapped_93.txt

coverm.0 <- read.table(file = "data/2023-03-15_coverM_on_drep96/output/coverM_reads_mapped_0.txt", sep = " ")
coverm.0 <- coverm.0[ ,c(7,9,14)]
colnames(coverm.0) <- c("sample","reads.mapped","reads.total")
coverm.0$perc.mapped <- coverm.0$reads.mapped / coverm.0$reads.total * 100
coverm.0$sample <- sub(",","",coverm.0$sample)
head(coverm.0)

coverm.93 <- read.table(file = "data/2023-03-15_coverM_on_drep96/output/coverM_reads_mapped_93.txt", sep = " ")
coverm.93 <- coverm.93[ ,c(7,9,14)]
colnames(coverm.93) <- c("sample","reads.mapped","reads.total")
coverm.93$perc.mapped <- coverm.93$reads.mapped / coverm.93$reads.total * 100
coverm.93$sample <- sub(",","",coverm.93$sample)
head(coverm.93)

coverm.reads.mapped <- list("ANI_0" = coverm.0, "ANI_93" = coverm.93)

saveRDS(object = coverm.reads.mapped, file = "data/2023-03-15_coverM_on_drep96/processed/reads_mapped.rds")

# ---- quick look ----

par(mfrow=c(2,1))
hist(coverm.0$perc.mapped)
hist(coverm.93$perc.mapped)
