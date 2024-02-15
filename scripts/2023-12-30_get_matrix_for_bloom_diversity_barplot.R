# RRR
# what proportion of all mags with seasonality had more/less diverse bloom patterns?
# barplot split by taxonomy

library(data.table)

seas <- fread(file = "data/2023-12-06_abundance_seasonality_analysis/genome_fft_stats_and_bloom_diversity.tsv.gz")
tax <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv")

# ---- format ----

seas <- seas[abund.is.seasonal == TRUE & nuc.div.is.seasonal == TRUE, .(genome,blooms.are)]

tax <- tax[ ,.(genome, domain, phylum, class, order, family, genus, species, mean.abund)]

seas <- merge(x = tax, y = seas, by = "genome", all.x = F, all.y = TRUE)

# ---- summarize by phylum ----

abundance.order <- seas[ , .(sum.mean.abund = sum(mean.abund)), by = phylum]
abundance.order <- abundance.order[order(sum.mean.abund, decreasing = TRUE)]
abundance.order[ ,bar.order := 1:nrow(abundance.order)]
abundance.order[ ,display.name := phylum]
abundance.order[sum.mean.abund < .5, display.name := "Other"]
abundance.order[sum.mean.abund < .5, bar.order := 8]

tot.genomes <- seas[ ,.(tot.genomes = .N), by = phylum]
tot.genomes <- merge(x = tot.genomes, y = abundance.order, by = "phylum")

tot.other.genomes <- sum(tot.genomes[display.name == "Other", tot.genomes])

phy <- seas[ , .(num.genomes = .N, mean.abund = sum(mean.abund, na.rm = TRUE)), by = .(phylum, blooms.are)]

phy <- merge(x = phy, y = tot.genomes, by = "phylum")

phy <- phy[ , .(num.genomes = sum(num.genomes), tot.genomes = sum(tot.genomes), mean.abund = sum(mean.abund)), by = .(display.name, blooms.are, bar.order)]

phy[display.name == "Other", tot.genomes := tot.other.genomes]

phy <- phy[ ,perc.genomes := num.genomes / tot.genomes * 100]

plotting.matrix <- dcast(data = phy, formula = blooms.are ~ display.name, value.var = "perc.genomes")

plotting.matrix <- as.matrix(plotting.matrix, rownames = T)

row.order <- c(3,2,1)
col.order <- merge(x = data.table("display.name" = colnames(plotting.matrix), "colnumber" = 1:ncol(plotting.matrix)), y = abundance.order[1:8,], by = "display.name")
col.order <- col.order[order(bar.order)]
col.order <- col.order$colnumber

plotting.matrix <- plotting.matrix[row.order, col.order]

# # check
# barplot(plotting.matrix, legend.text = TRUE)

saveRDS(object = plotting.matrix, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/bloom_diversity_summary_data.rds")
