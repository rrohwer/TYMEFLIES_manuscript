# RRR
# get the table of reads mapped plus season for making the plot
# fig 1c

library(data.table)

mapped <- readRDS("data/2023-03-15_coverM_on_drep96/processed/reads_mapped.rds")
sample.key <- fread("data/2023-12-08_combined_genome_info/sample_key.tsv")


mapped <- as.data.table(mapped)
mapped <- mapped[ ,.("sample" = ANI_93.sample, ANI_93.perc.mapped)]

sample.key <- sample.key[ ,.(sample, season)]

mapped <- merge(x = sample.key, y = mapped, by = "sample")

# boxplot wants a named list as input
mapped.list <- list()
for (s in unique(mapped$season)){
  my.list <- mapped[season == s, ANI_93.perc.mapped]
  my.list <- list(my.list)
  names(my.list) <- s
  mapped.list <- c(mapped.list, my.list)
}

# match order
names(mapped.list)
mapped.list <- c(mapped.list[6], mapped.list[1:5])

# # check
# boxplot(mapped.list, pch = NA)
# stripchart(x = mapped.list, vertical = T, method = "jitter", jitter = .25, add = T)

saveRDS(object = mapped.list, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_1_data/percent_mapped.rds")