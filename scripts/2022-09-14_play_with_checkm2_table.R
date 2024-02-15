x <- read.delim("data/2022-09-14_troubleshoot_checkm2/quality_report.tsv")

summary(x$Completeness)

hist(x$Completeness)

sum(x$Completeness > 50)

sum(x$Contamination < 10)

sum(x$Completeness > 50 & x$Contamination < 10)

