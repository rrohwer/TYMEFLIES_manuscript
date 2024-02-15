# RRR

library(data.table)
library(lubridate)
library(ggplot2)
library(patchwork)

# ---- do with GO-terms from interproscan ----

# fread and fwrite would be faster for large files, and can still zip IF Open MP enabled!
genes <- fread(file = "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39_gene_info.tsv")

tax <- readRDS(file = "data/2022-11-15_bin_stats/all_bins.rds")

anns <- fread(file = "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39-1Nlookup.tsv", sep = "\t", quote = "", fill = T)
colnames(anns) <- c("gene","MD5","length","analysis","signature.accession","signature.description","start","stop","e.value","match.status","run.date","interpro.accession","interpro.description","GO.annotation","pathways.annotation")
anns <- anns[GO.annotation != "" & GO.annotation != "-", .(gene, length, interpro.accession, interpro.description, GO.annotation)]

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
genes[ ,invasion := factor(invasion, levels = c("none","spiny","zebra"))]
colnames(genes)
colnames(anns)

length(unique(genes$gene)) # not unique because 471 observations of each gene
length(unique(anns$gene)) # not one line per gene
# multiple GO term lines for each annotation (pfam, tigerfam, etc) that had one
# multiple GO terms per line because there are different levels of GO?

anns <- anns[!duplicated(gene)] # can't merge with duplicate gene names, just take first one for now but how do I choose??

genes <- merge(x = genes, y = anns, by = "gene", all = T) 

genes$mcdonald.kreitman <- genes$pNpS_variants / genes$dNdS_substitutions # < 1 positive selection, >1 negative or balancing selection

hist(genes[ ,mcdonald.kreitman], breaks = 100)
hist(genes[ ,mcdonald.kreitman], breaks = 1000, xlim = c(0,10))
hist(log(genes[ ,mcdonald.kreitman]), breaks = 100)

genes$log.MK <- log(genes$mcdonald.kreitman)

# OK, how would you summarize by cogs under selection?
summary <- genes[log.MK > 1, .(.N), by = .(interpro.description)] # and add "COG" column to by
summary <- summary[order(-N)]

par(mar = c(20,3.5, 1, 0))
bar.spots <- barplot(height = summary$N[2:50], names.arg = rep("",49))
text(x = bar.spots, y = -2, labels = summary$interpro.description[2:50], xpd = NA, srt = 90, cex = .7, adj = 1)
mtext(text = "number genes", side = 2, line = 2.5)


# ---- do with COG terms from eggnog ----

# fread and fwrite would be faster for large files, and can still zip IF Open MP enabled!
genes <- fread(file = "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39_gene_info.tsv")

tax <- readRDS(file = "data/2022-11-15_bin_stats/all_bins.rds")

cogdefs <- fread(file = "data/2023-04-30_COG_info/COG_descriptions.csv")
colnames(cogdefs) <- c("COG_category", "COG_color", "COG_description")

anns <- fread(file = "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39.emapper.annotations", sep = "\t", quote = "", fill = T, skip = 4)
colnames(anns)[1] <- "query"
anns <- anns[-((nrow(anns) - 2):nrow(anns))] # stupid stats listed at end of file
unique(anns$COG_category)
length(unique(anns$COG_category)) # just 102 of them!
unique(anns$Description) # OK these are more detailed than the category
unique(anns$Preferred_name) # also just a handful
unique(anns$GOs) # ok of the ones that have it, there's a shit ton listed in a row
unique(anns$EC) # 427 total... these are just numbers of some sort
anns <- anns[COG_category != "-", .(query, COG_category, Description)]

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
genes[ ,invasion := factor(invasion, levels = c("none","spiny","zebra"))]
colnames(genes)
colnames(anns)

length(unique(genes$gene)) # not unique because 471 observations of each gene
length(unique(anns$query)) # ONE PER GENE AT LEAST

genes <- merge(x = genes, y = anns, by.x = "gene", by.y = "query", all = T) 

genes$mcdonald.kreitman <- genes$pNpS_variants / genes$dNdS_substitutions # < 1 positive selection, >1 negative or balancing selection

hist(genes[ ,mcdonald.kreitman], breaks = 100)
hist(genes[ ,mcdonald.kreitman], breaks = 1000, xlim = c(0,10))
hist(log(genes[ ,mcdonald.kreitman]), breaks = 100)

genes$log.MK <- log(genes$mcdonald.kreitman)

# OK, how would you summarize by cogs under selection?
target <- genes[log.MK > 1, .(.N), by = .(COG_category)] # and add "COG" column to by
target <- target[order(-N)]
target[is.na(COG_category), COG_category := "?"]

par(mar = c(4,3.5, 1, 0))
bar.spots <- barplot(height = target$N, names.arg = NULL)
text(x = bar.spots, y = -50, labels = target$COG_category, xpd = NA, srt = 90, cex = .5, adj = 1)
mtext(text = "number genes", side = 2, line = 2.5)
mtext(text = "COG category", side = 1, line = 2.5)
mtext(text = "Overall", side = 3, line = 0)

# OK, what to do about the duplicate (multi-digit) COG categories?
target
multis <- target[nchar(COG_category) > 1]
target <- target[nchar(COG_category) == 1]
# not sure how to add the multiple-category ones in.. remove for now

target <- merge(x = target, y = cogdefs, by = "COG_category", all.x = T, all.y = F)
target <- target[order(-N)]
color.key <- data.table("COG_description" = )


par(mar = c(4,3.5, 1, 0))
bar.spots <- barplot(height = target$N, names.arg = NULL)
text(x = bar.spots, y = -100, labels = target$COG_category, xpd = NA, cex = 1, adj = .5)
mtext(text = "Number Genes", side = 2, line = 2.5)
mtext(text = "COG category", side = 1, line = 2.5)
mtext(text = "Overall", side = 3, line = 0)

