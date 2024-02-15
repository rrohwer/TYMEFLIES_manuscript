# RRR

# is it microniches allowing multiple existing strains to emerge,
# or niche expansion occurring in real time?
# ecology or evolution?

# assumption with microniches- diversity is standing, but below detection limit
# this we can try to test!

library(data.table)

options(scipen = 9999)

genome <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")

abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")

genome <- merge(x = genome, y = abund, by.x = c("sample","genome"), by.y = c("Sample","Genome"), all = T)

genome <- genome[!is.na(nucl_diversity)]

plot(x = genome$nucl_diversity, y = genome$`Relative Abundance (%)`)
summary(genome$`Relative Abundance (%)`)

hist(genome$nucl_diversity, breaks = 100)
plot(x = genome[nucl_diversity < .01, nucl_diversity], y = genome[nucl_diversity < .01, `Relative Abundance (%)`])
plot(x = genome[nucl_diversity < .001, nucl_diversity], y = genome[nucl_diversity < .001, `Relative Abundance (%)`])

plot(x = genome[`Relative Abundance (%)` > 0 & !is.na(nucl_diversity), nucl_diversity], y = genome[`Relative Abundance (%)` > 0 & !is.na(nucl_diversity), `Relative Abundance (%)`])

plot(x = genome[ ,`Relative Abundance (%)`], y = genome[ ,nucl_diversity])
plot(x = genome[`Relative Abundance (%)` < 2, `Relative Abundance (%)`], y = genome[`Relative Abundance (%)` < 2, nucl_diversity])
plot(x = genome[`Relative Abundance (%)` < .01, `Relative Abundance (%)`], y = genome[`Relative Abundance (%)` < .01, nucl_diversity])

# it does not look like there's a "cutoff" where below some abundance you don't see the full nucleotide diversity
# this supports the niche expansion idea

# check: do pops have the same SNVs over time
# if niche expansion, would expect the SNVs to be different each season
# compare gives you popANI between 2 samples
# could also look at total SNV positions, do the SNVs in question change?
# compare generates a file pooled_SNV_info.tsv, info on each SNV across all samples

# So table with ALL SNVS in rows, all samples in columns, pres/abs for sample
# Over time, are SNVs lingering or changing?