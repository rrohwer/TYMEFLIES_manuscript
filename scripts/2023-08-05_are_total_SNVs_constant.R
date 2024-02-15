# RRR
# are total SNVs/sample somewhat constant (would suggest using Manhattan distances, like a max SNV load similar to a max species load)
# or are they variable (would suggest using Euclidean distances, since unlike comm ecol there is no unrealistically rich SNV space)

library(data.table)

genome.info <- fread(input = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
my.genome <- "ME2011-09-21_3300043464_group3_bin69"
genome.info <- genome.info[genome == my.genome]

genome.info[ ,above.breadth.cutoff := (breadth / breadth_expected) >= .5, ]

genome.info$divergent_site_count
hist(genome.info[above.breadth.cutoff == TRUE, divergent_site_count])
plot(genome.info$date, genome.info$divergent_site_count, type = "l", col = "red")
par(new = TRUE)
plot(genome.info$date, genome.info$coverage, type = "l", col = "blue")
# It's variable but it tracks the abundance very well

genome.info <- genome.info[above.breadth.cutoff == TRUE, sites.per.coverage := divergent_site_count / coverage]

plot(genome.info$date, genome.info$sites.per.coverage, type = "l")
# so although there is not an increase in NEW SNV sites, there is an increase in TOTAL SNV sites when normalized by coverage

# maybe I should normalize NEW sites by coverage too.... YES. but do that later

# but in any case, I think I need to use euclidean distances, the total SNVs is not constant
hist(genome.info$sites.per.coverage)
