# RRR
# forget singleM for now. 
# Instead, fill fig 1 with a comparison of rank abundance from 16S and MAGs
# Also can highlight acI abundance in it
# later can simply swap out the 16S side for singleM results


library(data.table)

library(limony)
data("limony")

genomes <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_sample.tsv.gz")

# ----

phylum.16S <- limony$av$Phylum
phylum.16S <- data.table("phylum" = rownames(phylum.16S), "Mean.Abund" = apply(X = phylum.16S, MARGIN = 1, FUN = mean))
phylum.16S <- phylum.16S[order(Mean.Abund, decreasing = T)]
save.other <- sum(phylum.16S[-c(1:7), Mean.Abund])
phylum.16S <- rbind(phylum.16S[1:7, ], data.table("phylum" = "Other", "Mean.Abund" = save.other))

order.16S <- limony$av$Order
order.16S <- data.table("order" = rownames(order.16S), "Mean.Abund" = apply(X = order.16S, MARGIN = 1, FUN = mean))
order.16S <- order.16S[order(Mean.Abund, decreasing = T)]
order.16S # Frankiales is most abundant order
nano.16S <- order.16S[order == "Frankiales", Mean.Abund]

phylum.mags <- genomes[lesser.duplicate == FALSE ,.(genome, phylum, sample, adj.abund.perc)]
phylum.mags <- phylum.mags[ , .(adj.abund.perc = sum(adj.abund.perc)), by = .(phylum, sample)]
phylum.mags <- phylum.mags[ , .(Mean.Abund = mean(adj.abund.perc)), by = .(phylum)]
phylum.mags <- phylum.mags[order(Mean.Abund, decreasing = T)]
save.other <- sum(phylum.mags[-c(1:7), Mean.Abund])
phylum.mags <- rbind(phylum.mags[1:7, ], data.table("phylum" = "Other", "Mean.Abund" = save.other))

order.mags <- genomes[lesser.duplicate == FALSE ,.(genome, order, sample, adj.abund.perc)]
order.mags <- order.mags[ , .(adj.abund.perc = sum(adj.abund.perc)), by = .(order, sample)]
order.mags <- order.mags[ , .(Mean.Abund = mean(adj.abund.perc)), by = .(order)]
order.mags <- order.mags[order(Mean.Abund, decreasing = T)]
order.mags # Nanopelagicales is most abundant order again
nano.mags <- order.mags[order == "Nanopelagicales", Mean.Abund]

# match names and put in 16S order
phylum.16S[phylum == "Chloroflexi", phylum := "Chloroflexota"]
phylum.mags <- phylum.mags[c(3,1,2,4,5,6,7,8)]
all.equal(phylum.16S$phylum, phylum.mags$phylum)

# add nano data for stacked bar
phylum.16S$nano.mean.abund <- c(0,nano.16S,0,0,0,0,0,0)
phylum.16S[phylum == "Actinobacteriota", Mean.Abund := Mean.Abund - nano.16S]
phylum.mags$nano.mean.abund <- c(0,nano.mags,0,0,0,0,0,0)
phylum.mags[phylum == "Actinobacteriota", Mean.Abund := Mean.Abund - nano.mags]

# save as matrixes for barplot input
phylum.16S <- as.matrix(x = phylum.16S, rownames = T)
phylum.mags <- as.matrix(x = phylum.mags, rownames = T)
phylum.16S <- t(phylum.16S)
phylum.mags <- t(phylum.mags)
phylum.16S <- phylum.16S[ ,ncol(phylum.16S):1]
phylum.mags <- phylum.mags[ ,ncol(phylum.mags):1]

# # test
# par(mar = c(3,8,1,1))
# barplot(phylum.16S, horiz = T, las = 2)
# barplot(phylum.mags, horiz = T, las = 2)
# all.equal(colnames(phylum.16S),colnames(phylum.mags))

# save ----

saveRDS(object = phylum.16S, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_1_data/16S_rank_abund.rds")
saveRDS(object = phylum.mags, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_1_data/MAG_rank_abund.rds")


# check Nanopelagicaceae mean abund as well, for text: ----

family.mags <- genomes[lesser.duplicate == FALSE ,.(genome, family, sample, adj.abund.perc)]
family.mags <- family.mags[ , .(adj.abund.perc = sum(adj.abund.perc)), by = .(family, sample)]
family.mags <- family.mags[ , .(Mean.Abund = mean(adj.abund.perc)), by = .(family)]
family.mags <- family.mags[order(Mean.Abund, decreasing = T)]
family.mags # Nanopelagicaceae 8.3428512070 Mean Abund
