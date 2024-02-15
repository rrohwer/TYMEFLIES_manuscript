# RRR
# what proportion of all mags had seasonality in pi and abundance?
# barplot split by taxonomy

library(data.table)

seas <- fread(file = "data/2023-12-06_abundance_seasonality_analysis/genome_fft_stats_and_bloom_diversity.tsv.gz")
tax <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv")

# ---- format ----

seas[ ,seasonality := "None"]
seas[abund.is.seasonal == TRUE & nuc.div.is.seasonal == TRUE, seasonality := "Both"]
seas[abund.is.seasonal == TRUE & nuc.div.is.seasonal == FALSE, seasonality := "Abundance"]
seas[abund.is.seasonal == FALSE & nuc.div.is.seasonal == FALSE, seasonality := "Diversity"]

seas <- seas[ ,.(genome, seasonality)]

tax <- tax[ ,.(genome, domain, phylum, class, order, family, genus, species, mean.abund)]

seas <- merge(x = tax, y = seas, by = "genome", all.x = F, all.y = TRUE)

length(unique(seas$genome)) # 1474 total genomes included in this analysis

# ---- summarize by phylum ----

abundance.order <- seas[ , .(sum.mean.abund = sum(mean.abund)), by = phylum]
abundance.order <- abundance.order[order(sum.mean.abund, decreasing = TRUE)]
abundance.order[ ,bar.order := 1:nrow(abundance.order)]
abundance.order[ ,display.name := phylum]
abundance.order[sum.mean.abund < 2, display.name := "Other"]
abundance.order[sum.mean.abund < 2, bar.order := 8]

tot.genomes <- seas[ ,.(tot.genomes = .N), by = phylum]
tot.genomes <- merge(x = tot.genomes, y = abundance.order, by = "phylum")

tot.other.genomes <- sum(tot.genomes[display.name == "Other", tot.genomes])

phy <- seas[ , .(num.genomes = .N, mean.abund = sum(mean.abund, na.rm = TRUE)), by = .(phylum, seasonality)]

phy <- merge(x = phy, y = tot.genomes, by = "phylum")

phy <- phy[ , .(num.genomes = sum(num.genomes), tot.genomes = sum(tot.genomes), mean.abund = sum(mean.abund)), by = .(display.name, seasonality, bar.order)]

phy[display.name == "Other", tot.genomes := tot.other.genomes]

phy <- phy[ ,perc.genomes := num.genomes / tot.genomes * 100]

plotting.matrix <- dcast(data = phy, formula = seasonality ~ display.name, value.var = "perc.genomes")

plotting.matrix <- as.matrix(plotting.matrix, rownames = T)

row.order <- c(4,1,2,3)
col.order <- merge(x = data.table("display.name" = colnames(plotting.matrix), "colnumber" = 1:ncol(plotting.matrix)), y = abundance.order[1:8,], by = "display.name")
col.order <- col.order[order(bar.order)]
col.order <- col.order$colnumber

plotting.matrix <- plotting.matrix[row.order, col.order]

# # check
# barplot(plotting.matrix, legend.text = TRUE)

saveRDS(object = plotting.matrix, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/seasonality_summary_data.rds")
