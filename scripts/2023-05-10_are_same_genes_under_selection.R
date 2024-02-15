# RRR

# Are the same genes repeatedly under selection, or are different genes w/in the cog category under selection

# ---- set up ----

library(data.table)
library(lubridate)
library(ggplot2)

genes <-"data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39_gene_info.tsv"
cogdefs <- "data/2023-04-30_COG_info/COG_descriptions.csv"
anns <- "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39.emapper.annotations"
genomes <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"
threads <- 1

genes <- fread(file = genes, nThread = threads)
cogdefs <- fread(file = cogdefs, nThread = threads)
anns <- fread(file = anns, sep = "\t", quote = "", fill = T, skip = 4, nThread = threads)
genomes <- readRDS(genomes)

# ---- merge data ----

colnames(anns)[1] <- "query"
anns <- anns[-((nrow(anns) - 2):nrow(anns))] # stupid stats listed at end of file
anns <- anns[COG_category != "-", .(query, COG_category, Description)]

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
genes[ ,invasion := factor(invasion, levels = c("none","spiny","zebra"))]

# colnames(genes)
# colnames(anns)
# 
# length(unique(genes$gene)) # not unique because 471 observations of each gene
# length(unique(anns$query)) # ONE PER GENE AT LEAST

genes <- merge(x = genes, y = anns, by.x = "gene", by.y = "query", all = T) 

# colnames(genes)
# colnames(cogdefs)

genes <- merge(x = genes, y = cogdefs, by = "COG_category", all.x = T, all.y = F)

genes <- merge(x = genes, y = genomes, by.x = "genome", by.y = "bin.full.name", all.x = T, all.y = F)

tax.label <- paste(genes[1, .(phylum, class, order, family, genus, species)], collapse = " ")

my.genome <- unique(genes$genome)

plot.label <- paste0(tax.label, " ", my.genome,".pdf")


# ---- ID selection ----

genes$mcdonald.kreitman <- genes$pNpS_variants / genes$dNdS_substitutions # < 1 positive selection, >1 negative or balancing selection

# hist(genes[ ,mcdonald.kreitman], breaks = 100)
# hist(genes[ ,mcdonald.kreitman], breaks = 1000, xlim = c(0,10))
# hist(log(genes[ ,mcdonald.kreitman]), breaks = 100)

genes$log.MK <- log(genes$mcdonald.kreitman)

# ---- times/gene ----

genes[is.na(COG_category), COG_category := "?"]
genes <- genes[nchar(COG_category) == 1]
length(unique(genes$gene))
nrow(unique(genes[ ,c("gene","COG_category")]))

sel.per.gene <- genes[abs(log.MK) >= 1.5 & is.finite(log.MK), .N, by = .(genome, gene, COG_category)]

hist(sel.per.gene$N, breaks = 100)
sel.per.gene <- sel.per.gene[order(-N)]
barplot(height = sel.per.gene[ ,N])
barplot(height = sel.per.gene[1:100, N])
boxplot(sel.per.gene[ ,N])

sel.per.gene.per.cog <- dcast(data = sel.per.gene, formula = gene ~ COG_category, value.var = "N")
boxplot(x = sel.per.gene.per.cog[ ,-1], lty = 1)

boxplot(x = sel.per.gene.per.cog[ ,-1], lty = 1, range = 0)
stripchart(x = sel.per.gene.per.cog[ ,-1], method = "jitter", vertical = T, add = T, pch = 19, col = adjustcolor("black",.2), jitter = c(.3,0))

stripchart(x = sel.per.gene.per.cog[ ,-1], method = "jitter", vertical = T, add = F, pch = 19, col = adjustcolor("black",.2), jitter = c(.3,0))

# So select genes are repeatedly under selection
# and in some categories the same genes are esp repeatedly under selection
# so could filter by how many times the gene was observed under selection too.

# ---- new genes per obs ----

sel <- genes[abs(log.MK) >= 1.5  & is.finite(log.MK), .(date, gene, COG_category, log.MK)]
sel <- sel[order(date)]
sel[ ,first.obs := !duplicated(gene)]
sel[ ,new.obs := sum(first.obs), by = date]
collect.overall <- sel[!duplicated(date), .(date,new.obs)]
collect.overall[ ,cum := cumsum(new.obs)]

ggplot(data = collect.overall, aes(x = date, y = cum))+
  geom_point()+
  ylab("New genes under selection")
  
sel <- genes[abs(log.MK) >= 1.5  & is.finite(log.MK), .(date, gene, COG_category, log.MK)]
sel <- sel[order(date)]
sel[ ,first.obs := !duplicated(gene), by = COG_category]
sel[ ,new.obs := sum(first.obs), by = .(date,COG_category)]
collect.cog <- sel[ , .N, by = .(date, COG_category, new.obs)]
collect.cog <- collect.cog[ ,-"N"]
collect.cog[ ,cum := cumsum(new.obs), by = COG_category]

ggplot(data = collect.cog, aes(x = date, y = cum))+
  geom_point()+
  ylab("New genes under selection")

ggplot(data = collect.cog[COG_category != "?" & COG_category != "S"], aes(x = date, y = cum))+
  geom_point()+
  ylab("New genes under selection")+
  facet_wrap(~COG_category)

# but accumulating new genes over time pretty consistently. except O... maybe all O's were in 1 sample, cool.


