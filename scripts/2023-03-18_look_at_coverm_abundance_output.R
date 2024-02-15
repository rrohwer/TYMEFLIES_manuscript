# RRR
# get the actual coverM output files into RDS format

library(data.table)

ani0 <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI0.txt")
ani93 <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt")

colnames(ani0)
head(ani0)

x <- readRDS(file = "data/2023-03-15_coverM_on_drep96/processed/reads_mapped.rds")
x <- x$ANI_0
y = ani0[Sample == "ME2013-02-02s1D0_3300042325" & Genome == "unmapped","Relative Abundance (%)"]
z = x[x$sample == "ME2013-02-02s1D0_3300042325", "perc.mapped"]
y + z # 100

# so the "unmapped" contigs are the reads that didn't map to anything. doh, could have gotten percent coverage from this file!


