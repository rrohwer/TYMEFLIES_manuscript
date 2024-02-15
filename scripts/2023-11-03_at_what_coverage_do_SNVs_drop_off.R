# RRR
# at what coverage do the number of SNVs called drop off?
# i.e. at what coverage should we exclude dates from the SNV analysis, because SNVs may be being called unreliably?

library(data.table)

genome.info <- fread(input = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")

genome.info[ ,above.breadth.cutoff := (breadth / breadth_expected) >= .5, ]

plot(genome.info$coverage_median, genome.info$divergent_site_count, type = "n")
points(genome.info$coverage_median, genome.info$divergent_site_count, col = adjustcolor("black",.1))

plot(genome.info[above.breadth.cutoff == TRUE, coverage_median], genome.info[above.breadth.cutoff == TRUE, divergent_site_count], type = "n")
points(genome.info[above.breadth.cutoff == TRUE, coverage_median], genome.info[above.breadth.cutoff == TRUE, divergent_site_count], col = adjustcolor("black",.1))

plot(genome.info[above.breadth.cutoff == TRUE & coverage_median < 200, coverage_median], genome.info[above.breadth.cutoff == TRUE & coverage_median < 200, divergent_site_count], type = "n")
points(genome.info[above.breadth.cutoff == TRUE & coverage_median < 200, coverage_median], genome.info[above.breadth.cutoff == TRUE & coverage_median < 200, divergent_site_count], col = adjustcolor("black",.4))

plot(genome.info[above.breadth.cutoff == TRUE & coverage_median < 50, coverage_median], genome.info[above.breadth.cutoff == TRUE & coverage_median < 50, divergent_site_count], type = "n")
points(genome.info[above.breadth.cutoff == TRUE & coverage_median < 50, coverage_median], genome.info[above.breadth.cutoff == TRUE & coverage_median < 50, divergent_site_count], col = adjustcolor("black",.4))
abline(v = 10)
abline(v = 5)

# median coverage == 10 seems like a good cutoff to remove those few high outliers
# but median coverage > 5 seems like a good cutoff to chop it before the bulk of the genomes decrease in SNVs
# maybe let's add in rel abund info?
# bumping breadth/breadth expected to .75 doesn't change the story that much.

abund.info <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
colnames(abund.info)[3] <- "rel.abund"

genome.info <- merge(x = abund.info, y = genome.info, by.x = c("Sample","Genome"), by.y = c("sample","genome"), all = T)

plot(genome.info[above.breadth.cutoff == TRUE & rel.abund > 0, coverage_median], genome.info[above.breadth.cutoff == TRUE & rel.abund > 0, divergent_site_count], type = "n")
points(genome.info[above.breadth.cutoff == TRUE & rel.abund > 0, coverage_median], genome.info[above.breadth.cutoff == TRUE & rel.abund > 0, divergent_site_count], col = adjustcolor("black",.1))

plot(genome.info[above.breadth.cutoff == TRUE & coverage_median < 50 & rel.abund > .01, coverage_median], genome.info[above.breadth.cutoff == TRUE & coverage_median < 50 & rel.abund > .01, divergent_site_count], type = "n")
points(genome.info[above.breadth.cutoff == TRUE & coverage_median < 50 & rel.abund > .01, coverage_median], genome.info[above.breadth.cutoff == TRUE & coverage_median < 50 & rel.abund > .01, divergent_site_count], col = adjustcolor("black",.4))
abline(v = 10)
abline(v = 5)
# still looks the same

plot(genome.info[above.breadth.cutoff == TRUE & coverage_median < 50 & rel.abund > .05, coverage_median], genome.info[above.breadth.cutoff == TRUE & coverage_median < 50 & rel.abund > .05, divergent_site_count], type = "n")
points(genome.info[above.breadth.cutoff == TRUE & coverage_median < 50 & rel.abund > .05, coverage_median], genome.info[above.breadth.cutoff == TRUE & coverage_median < 50 & rel.abund > .05, divergent_site_count], col = adjustcolor("black",.4))
abline(v = 10)
abline(v = 5)
# still looks the same, so abund is not really helping much.

library(ggplot2)
x <- genome.info[above.breadth.cutoff == TRUE & coverage_median < 50, .(coverage_median, divergent_site_count)]
ggplot(data = x, aes(y = divergent_site_count, x = coverage_median))+
  geom_boxplot(aes(group = coverage_median))

ggplot(data = x, aes(y = divergent_site_count, x = coverage_median))+
  geom_boxplot(aes(group = coverage_median))+
  coord_cartesian(ylim = c(0,250000))+
  geom_line(aes(x = 10.5), color = "red")

# AND how much are these cutoffs reducing the number of samples included?

nrow(genome.info) # 1,345,176 = 2855 genomes x 471 samples
sum(genome.info$rel.abund > 0) # 469,215 genomes where it's present at all by rel abund
sum((genome.info[!is.na(breadth) & !is.na(breadth_expected), breadth] / genome.info[!is.na(breadth) & !is.na(breadth_expected), breadth_expected]) >= .5) # 240,278 where it's above .5 breadth expected
sum((genome.info[!is.na(breadth) & !is.na(breadth_expected), breadth] / genome.info[!is.na(breadth) & !is.na(breadth_expected), breadth_expected]) >= .75) # 208,470 where it's above .75 breadth expected
nrow(genome.info[coverage_median > 5]) # 68,852 where coverage is above 5
nrow(genome.info[coverage_median > 10]) # 39,954 where coverage is above 10

# OK, so median genome coverage > 5 is a pretty harsh cutoff, and > 10 you throw out almost half the remaining sample to remove like 10 outlier high SNV samples
# No now that I look with boxplots, the distribution is more clear. And it is definitely in the slope of the drop-off still at 5
# I think inStrain says to trust SNVs if coverage is > 5, but that is the SNV coverage and we have to filter based on average genome coverage

# How much does this decrease the number of samples for our example genomes?

acI.B <- genome.info[Genome == "ME2011-09-21_3300043464_group3_bin69"]
nrow(acI.B)
nrow(acI.B[rel.abund > 0]) # 471
nrow(acI.B[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5]) # 471
nrow(acI.B[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5 & coverage_median > 5]) # 417 , so lost 54 sample dates
nrow(acI.B[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5 & coverage_median > 10]) # 307 , so lost 164 samples damn
acI.A <- genome.info[Genome == "ME2011-09-04_3300044729_group3_bin142"]
nrow(acI.A)
nrow(acI.A[rel.abund > 0]) # 470
nrow(acI.A[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5]) # 459
nrow(acI.A[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5 & coverage_median > 5]) # 301 , so lost 158 sample dates
nrow(acI.A[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5 & coverage_median > 10]) # 211 , so lost 248 samples damn
acI.C <- genome.info[Genome == "ME2016-07-20_3300033996_group7_bin32"]
nrow(acI.C)
nrow(acI.C[rel.abund > 0]) # 446
nrow(acI.C[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5]) # 363
nrow(acI.C[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5 & coverage_median > 5]) # 231 , so lost 132 sample dates
nrow(acI.C[!is.na(breadth) & !is.na(breadth_expected) & (breadth / breadth_expected) >= .5 & coverage_median > 10]) # 205 , so lost 158 samples damn

# final decision: include sample in SNV analysis when the genome is >= .5 breadth/breadth expected (as I've been doing) and when median coverage is > 10 