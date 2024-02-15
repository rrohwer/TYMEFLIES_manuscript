# RRR

# ---- set up ----

library(data.table)
library(lubridate)

genes <- fread(file = "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39_gene_info.tsv")

tax <- readRDS(file = "data/2022-11-15_bin_stats/all_bins.rds")

cogdefs <- fread(file = "data/2023-04-30_COG_info/COG_descriptions.csv")

anns <- fread(file = "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39.emapper.annotations", sep = "\t", quote = "", fill = T, skip = 4)

genomes <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

# ---- merge data ----

colnames(anns)[1] <- "query"
anns <- anns[-((nrow(anns) - 2):nrow(anns))] # stupid stats listed at end of file
anns <- anns[COG_category != "-", .(query, COG_category, Description)]

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
genes[ ,invasion := factor(invasion, levels = c("none","spiny","zebra"))]

colnames(genes)
colnames(anns)

length(unique(genes$gene)) # not unique because 471 observations of each gene
length(unique(anns$query)) # ONE PER GENE AT LEAST

genes <- merge(x = genes, y = anns, by.x = "gene", by.y = "query", all = T) 

colnames(genes)
colnames(cogdefs)

genes <- merge(x = genes, y = cogdefs, by = "COG_category", all.x = T, all.y = F)

genes <- merge(x = genes, y = genomes, by.x = "genome", by.y = "bin.full.name", all.x = T, all.y = F)

# ---- ID selection ----

genes$mcdonald.kreitman <- genes$pNpS_variants / genes$dNdS_substitutions # < 1 positive selection, >1 negative or balancing selection

hist(genes[ ,mcdonald.kreitman], breaks = 100)
hist(genes[ ,mcdonald.kreitman], breaks = 1000, xlim = c(0,10))
hist(log(genes[ ,mcdonald.kreitman]), breaks = 100)

genes$log.MK <- log(genes$mcdonald.kreitman)

sum(!(genes$SNS_N_count < 5 | genes$SNS_S_count < 5 | genes$SNV_N_count < 5 | genes$SNV_S_count < 5), na.rm = T) # **** should I do chi-square test instead of fishers exact test for these ones?****

get.fisher.p.value <- function(X){
  X <- matrix(X, nrow = 2, byrow = T)
  if(any(is.na(X))){
    X <- NA
  }else{
    X <- fisher.test(X)$p.value 
  }
  return(X)
}

add.fisher.p.value <- function(my.genes){
  my.mat <- as.matrix(my.genes[ ,.(SNV_N_count, SNV_S_count, SNS_N_count, SNS_S_count)])
  my.genes$MK.p.val <- apply(X = my.mat, MARGIN = 1, FUN = get.fisher.p.value)
  return(my.genes)
}

genes <- add.fisher.p.value(my.genes = genes)

length(which(genes$MK.p.val < .05))

# ---- overall selected COGs ----

pos <- genes[log.MK < -1.5 & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description,Broad_category)]
pos[is.na(COG_category), COG_category := "?"]

neg <- genes[log.MK > 1.5 & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description, Broad_category)]
neg[is.na(COG_category), COG_category := "?"]

# OK, what to do about the duplicate (multi-digit) COG categories?
# not sure how to add the multiple-category ones in.. remove for now *************
pos <- pos[nchar(COG_category) == 1]
neg <- neg[nchar(COG_category) == 1]

# And not sure whether to include the unknown functions *******
# ALSO, why are all the NA's not in the S category?
pos <- pos[COG_category != "?" & COG_category != "S"]
neg <- neg[COG_category != "?" & COG_category != "S"]

sel <- merge(x = pos, y = neg, by = c("COG_category","COG_Description","Broad_category"), all = T, suffixes = c("pos","neg"))

# ---- make overall plots ----

tax.label <- paste(genes[1, .(phylum, class, order, family, genus, species)], collapse = " ")

# ---- plot of all COG categories ----

col.key <- data.table("Broad_category" = c(NA,"Poorly Characterized","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                      "color" = c("grey70","grey70","skyblue3","seagreen","rosybrown3"))
sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)

index <- order(sel$Npos, decreasing = T)
sel <- sel[index]

bar.spots <- barplot(height = sel$Npos, names.arg = NULL, col = pos$color, axes = F, plot = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(sel$Nneg, na.rm = T), max(sel$Npos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot1.pdf", width = 6.5, height = 3)
par(mar = c(1.5,3.5, 3, 13))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = sel$Npos, names.arg = NULL, col = sel$color, axes = F, add = T)
barplot(height = -sel$Nneg, names.arg = NULL, col = sel$color, axes = F, add = T)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 15, labels = sel$COG_category, xpd = NA, srt = 0, adj = .5, cex = .7)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes (count)", side = 2, line = 2.75, cex = .7)
mtext(text = "Positive Selection", side = 2, line = 1.75, at = y.lim[2] / 2, cex = .7)
mtext(text = "Purifying Selection", side = 2, line = 1.75, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .4, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1.25, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2.15, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 20, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:length(bar.spots)), labels = paste(sel$COG_category,"=",sel$COG_Description), 
     xpd = NA, adj = 0, cex = .5)
mtext(text = col.key[-(1:2),Broad_category], col = col.key[-(1:2),color],
      side = 3,line = c(1.2,0.5,-.1), cex = .7, at = x.lim[2] + x.lim[2] / 20, adj = 0)
dev.off()

# ---- plot of all metabolism COGs ----

sel <- sel[Broad_category == "Metabolism"]

metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                             "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                             "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                           "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                           "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                       "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, and catabolism",
                                       "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
sel <- merge(x = sel, y = metab.colors, by = "COG_Description")

index <- order(sel$Npos, decreasing = T)
sel <- sel[index]

bar.spots <- barplot(height = sel$Npos, names.arg = NULL, col = sel$color, axes = F, plot = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(sel$Nneg, na.rm = T), max(sel$Npos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot2.pdf", width = 6.5, height = 3)
par(mar = c(1.5,4.5, 2, 9.75))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = sel$Npos, names.arg = NULL, col = sel$metab.color, axes = F, add = T)
barplot(height = -sel$Nneg, names.arg = NULL, col = sel$metab.color, axes = F, add = T)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 18, labels = sel$COG_category, xpd = NA, srt = 0, adj = .5, cex = .7)
ticks <- axis(side = 2, line = .5, labels = F, tck = -.02)
axis(side = 2, line = 0, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes (count)", side = 2, line = 3.75, cex = .7)
mtext(text = "Positive Selection", side = 2, line = 2.5, at = y.lim[2] / 2, cex = .7)
mtext(text = "Purifying Selection", side = 2, line = 2.5, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .4, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = .25, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 1.25, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 10, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:length(bar.spots)), labels = paste(sel$COG_category,"=",sel$label), 
     xpd = NA, adj = 0, cex = .6, col = sel$metab.color)
dev.off()

# ---- breakdown by season ----

pos <- genes[log.MK < -1.5  & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description, Broad_category, season)]
pos[is.na(COG_category), COG_category := "?"]

neg <- genes[log.MK > 1.5 & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description, Broad_category, season)]
neg[is.na(COG_category), COG_category := "?"]

# OK, what to do about the duplicate (multi-digit) COG categories?
# not sure how to add the multiple-category ones in.. remove for now *************
pos <- pos[nchar(COG_category) == 1]
neg <- neg[nchar(COG_category) == 1]

# And not sure whether to include the unknown functions *******
# ALSO, why are all the NA's not in the S category?
pos <- pos[COG_category != "?" & COG_category != "S"]
neg <- neg[COG_category != "?" & COG_category != "S"]

sel <- merge(x = pos, y = neg, by = c("COG_category","COG_Description","Broad_category","season"), all = T, suffixes = c("pos","neg"))

# ---- plot all COG categories seasonally ----

col.key <- data.table("Broad_category" = c(NA,"Poorly Characterized","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                      "color" = c("grey70","grey70","skyblue3","seagreen","rosybrown3"))
sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)

index <- order(sel$Npos, decreasing = T)
sel <- sel[index]

pos <- dcast(data = sel, formula = COG_category + color + COG_Description ~ season, value.var = "Npos")
key <- pos[ ,1:3]
pos <- as.matrix(pos[ ,-(2:3)], rownames = T)

neg <- dcast(data = sel, formula = COG_category ~ season, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)

all.equal(row.names(pos), row.names(neg))
all.equal(colnames(pos), colnames(neg))
all.equal(row.names(pos), key$COG_category)

bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, axes = F, plot = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(neg, na.rm = T), max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot3.pdf", width = 6.5, height = 3)
par(mar = c(1.25,3.75, 2.75, 9))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, space = c(0,4), ann = F, axisnames = F, col = key$color, add = T)
barplot(height = -neg, beside = T, names.arg = NULL, axes = F, space = c(0,4), ann = F, axisnames = F, col = key$color, add = T)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .3)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes (count)", side = 2, line = 3, cex = .7, at = 0)
mtext(text = "Positive Selection", side = 2, line = 2, at = 0, cex = .7, adj = 0)
mtext(text = "Purifying Selection", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1.25, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 25, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$COG_Description), 
     xpd = NA, adj = 0, cex = .3)
mtext(text = col.key[-(1:2),Broad_category], col = col.key[-(1:2),color],
      side = 3,line = c(1.2,0.5,-.1), cex = .7, at = x.lim[2] + x.lim[2] / 25, adj = 0)
dev.off()


# ---- plot all metabolism COGs seasonally ----

sel <- sel[Broad_category == "Metabolism"]

metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                                                 "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                                                 "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                           "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                           "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                       "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, catabolism",
                                       "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
sel <- merge(x = sel, y = metab.colors, by = "COG_Description")

pos <- dcast(data = sel, formula = COG_category + metab.color + COG_Description + label ~ season, value.var = "Npos")
key <- pos[ ,1:4]
pos <- as.matrix(pos[ ,-(2:4)], rownames = T)

neg <- dcast(data = sel, formula = COG_category ~ season, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)

all.equal(row.names(pos), row.names(neg))
all.equal(colnames(pos), colnames(neg))
all.equal(row.names(pos), key$COG_category)


bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(neg, na.rm = T), max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot4.pdf", width = 6.5, height = 3)
par(mar = c(1.25,3.75, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
barplot(height = -neg, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .6)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes (count)", side = 2, line = 3, cex = .7, at = 0)
mtext(text = "Positive Selection", side = 2, line = 2, at = y.lim[2] / 2, cex = .7, adj = .5)
mtext(text = "Purifying Selection", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)
dev.off()


# ---- plot all metabolism COGs seasonally- pos stacked ----

bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F)
x.lim <- c(0, max(bar.spots) + bar.spots[1])
y.lim <- c(0, max(colSums(pos, na.rm = T), na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot5.pdf", width = 6.5, height = 3)
par(mar = c(1.25,1.75, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
text(x = bar.spots, y = 0 - y.lim[2] / 10, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = -1, labels = F, tck = -.02)
axis(side = 2, line = -1.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Positive Selection (count)", side = 2, line = .75, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 50, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$label), 
     xpd = NA, adj = 0, cex = .7, col = key$metab.color)
dev.off()

# ---- plot all metabolism COGs seasonally- pos beside ----

bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(0, max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot6.pdf", width = 6.5, height = 3)
par(mar = c(1.25,2, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
text(x = bar.spots, y = 0 - y.lim[2] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .6)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Positive Selection (count)", side = 2, line = 1.25, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)
dev.off()





























# ---- breakdown by year ----

pos <- genes[log.MK < -1.5  & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description, Broad_category, year)]
pos[is.na(COG_category), COG_category := "?"]

neg <- genes[log.MK > 1.5 & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description, Broad_category, year)]
neg[is.na(COG_category), COG_category := "?"]

# OK, what to do about the duplicate (multi-digit) COG categories?
# not sure how to add the multiple-category ones in.. remove for now *************
pos <- pos[nchar(COG_category) == 1]
neg <- neg[nchar(COG_category) == 1]

# And not sure whether to include the unknown functions *******
# ALSO, why are all the NA's not in the S category?
pos <- pos[COG_category != "?" & COG_category != "S"]
neg <- neg[COG_category != "?" & COG_category != "S"]

sel <- merge(x = pos, y = neg, by = c("COG_category","COG_Description","Broad_category","year"), all = T, suffixes = c("pos","neg"))

# ---- plot all COG categories annually ----

col.key <- data.table("Broad_category" = c(NA,"Poorly Characterized","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                      "color" = c("grey70","grey70","skyblue3","seagreen","rosybrown3"))
sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)

index <- order(sel$Npos, decreasing = T)
sel <- sel[index]

pos <- dcast(data = sel, formula = COG_category + color + COG_Description ~ year, value.var = "Npos")
key <- pos[ ,1:3]
pos <- as.matrix(pos[ ,-(2:3)], rownames = T)

neg <- dcast(data = sel, formula = COG_category ~ year, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)

all.equal(row.names(pos), row.names(neg))
all.equal(colnames(pos), colnames(neg))
all.equal(row.names(pos), key$COG_category)

bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, axes = F, plot = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(neg, na.rm = T), max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot13.pdf", width = 12.5, height = 3)
par(mar = c(1.25,3.75, 2.75, 9))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, space = c(0,4), ann = F, axisnames = F, col = key$color, add = T, border = NA)
barplot(height = -neg, beside = T, names.arg = NULL, axes = F, space = c(0,4), ann = F, axisnames = F, col = key$color, add = T, border = NA)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .3)
# text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes (count)", side = 2, line = 3, cex = .7, at = 0)
mtext(text = "Positive Selection", side = 2, line = 2, at = 0, cex = .7, adj = 0)
mtext(text = "Purifying Selection", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1.25, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 25, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$COG_Description), 
     xpd = NA, adj = 0, cex = .3)
mtext(text = col.key[-(1:2),Broad_category], col = col.key[-(1:2),color],
      side = 3,line = c(1.2,0.5,-.1), cex = .7, at = x.lim[2] + x.lim[2] / 25, adj = 0)
mtext(text = colnames(pos), side = 3, line = -.2, at = bar.spots[floor(nrow(bar.spots)/2), ], cex = .5)
dev.off()


# ---- plot all metabolism COGs seasonally ----

sel <- sel[Broad_category == "Metabolism"]

metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                                                 "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                                                 "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                           "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                           "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                       "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, catabolism",
                                       "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
sel <- merge(x = sel, y = metab.colors, by = "COG_Description")

pos <- dcast(data = sel, formula = COG_category + metab.color + COG_Description + label ~ year, value.var = "Npos")
key <- pos[ ,1:4]
pos <- as.matrix(pos[ ,-(2:4)], rownames = T)

neg <- dcast(data = sel, formula = COG_category ~ year, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)

all.equal(row.names(pos), row.names(neg))
all.equal(colnames(pos), colnames(neg))
all.equal(row.names(pos), key$COG_category)


bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(neg, na.rm = T), max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot14.pdf", width = 12.5, height = 3)
par(mar = c(1.25,3.75, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
barplot(height = -neg, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .6)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes (count)", side = 2, line = 3, cex = .7, at = 0)
mtext(text = "Positive Selection", side = 2, line = 2, at = y.lim[2] / 2, cex = .7, adj = .5)
mtext(text = "Purifying Selection", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)
mtext(text = colnames(pos), side = 3, line = -.2, at = bar.spots[floor(nrow(bar.spots)/2), ], cex = .5)
dev.off()


# ---- plot all metabolism COGs annually- pos stacked ----

bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F)
x.lim <- c(0, max(bar.spots) + bar.spots[1])
y.lim <- c(0, max(colSums(pos, na.rm = T), na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot15.pdf", width = 6.5, height = 3)
par(mar = c(1.25,1.75, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
# text(x = bar.spots, y = 0 - y.lim[2] / 10, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = -1, labels = F, tck = -.02)
axis(side = 2, line = -1.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Positive Selection (count)", side = 2, line = .75, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 50, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$label), 
     xpd = NA, adj = 0, cex = .7, col = key$metab.color)
mtext(text = colnames(pos), side = 3, line = -.2, at = bar.spots, cex = .5)

# *********WHY IS THIS NOT WORKING????********

dev.off()

# ****** left off here adding seasonally ones *********

# ---- plot all metabolism COGs annually- pos beside ----

bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(0, max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot6.pdf", width = 6.5, height = 3)
par(mar = c(1.25,2, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
text(x = bar.spots, y = 0 - y.lim[2] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .6)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Positive Selection (count)", side = 2, line = 1.25, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)
dev.off()
# ---- Now repeat, but by percent of the COG category ----

# ---- overall selected COGs ----

pos <- genes[log.MK < -1.5 & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description,Broad_category)]
pos[is.na(COG_category), COG_category := "?"]

neg <- genes[log.MK > 1.5 & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description, Broad_category)]
neg[is.na(COG_category), COG_category := "?"]

# OK, what to do about the duplicate (multi-digit) COG categories?
# not sure how to add the multiple-category ones in.. remove for now *************
pos <- pos[nchar(COG_category) == 1]
neg <- neg[nchar(COG_category) == 1]

# And not sure whether to include the unknown functions *******
# ALSO, why are all the NA's not in the S category?
pos <- pos[COG_category != "?" & COG_category != "S"]
neg <- neg[COG_category != "?" & COG_category != "S"]

sel <- merge(x = pos, y = neg, by = c("COG_category","COG_Description","Broad_category"), all = T, suffixes = c("pos","neg"))

# tot.per.COG <- anns[ ,.(total.in.COG = .N), by = COG_category] # this is total COGs present, but for combining seasons seasons should be total observations of COGs
tot.per.COG <- genes[ ,.(total.in.COG = .N), by = COG_category] # consider adjusting per-season? ***********************

# Make it a percent of the total annotations in that COG category, since some COGs may be more abundant in the genomes
sel <- merge(x = sel, y = tot.per.COG, by = "COG_category", all.x = T, all.y = F)
sel[ ,`:=`(Npos = Npos / total.in.COG * 100, Nneg = Nneg / total.in.COG * 100)]

# ---- make overall plots ----

tax.label <- paste(genes[1, .(phylum, class, order, family, genus, species)], collapse = " ")

# ---- plot of all COG categories ----

col.key <- data.table("Broad_category" = c(NA,"Poorly Characterized","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                      "color" = c("grey70","grey70","skyblue3","seagreen","rosybrown3"))
sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)

index <- order(sel$Npos, decreasing = T)
sel <- sel[index]

bar.spots <- barplot(height = sel$Npos, names.arg = NULL, col = pos$color, axes = F, plot = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(sel$Nneg, na.rm = T), max(sel$Npos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot7.pdf", width = 6.5, height = 3)
par(mar = c(1.5,3.5, 3, 13))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = sel$Npos, names.arg = NULL, col = sel$color, axes = F, add = T)
barplot(height = -sel$Nneg, names.arg = NULL, col = sel$color, axes = F, add = T)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 15, labels = sel$COG_category, xpd = NA, srt = 0, adj = .5, cex = .7)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes per COG category (%)", side = 2, line = 2.75, cex = .7)
mtext(text = "Positive Selection", side = 2, line = 1.75, at = y.lim[2] / 2, cex = .7)
mtext(text = "Purifying Selection", side = 2, line = 1.75, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .4, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1.25, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2.15, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 20, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:length(bar.spots)), labels = paste(sel$COG_category,"=",sel$COG_Description), 
     xpd = NA, adj = 0, cex = .5)
mtext(text = col.key[-(1:2),Broad_category], col = col.key[-(1:2),color],
      side = 3,line = c(1.2,0.5,-.1), cex = .7, at = x.lim[2] + x.lim[2] / 20, adj = 0)
dev.off()

# ---- plot of all metabolism COGs ----

sel <- sel[Broad_category == "Metabolism"]

metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                                                 "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                                                 "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                           "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                           "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                       "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, and catabolism",
                                       "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
sel <- merge(x = sel, y = metab.colors, by = "COG_Description")

index <- order(sel$Npos, decreasing = T)
sel <- sel[index]

bar.spots <- barplot(height = sel$Npos, names.arg = NULL, col = sel$color, axes = F, plot = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(sel$Nneg, na.rm = T), max(sel$Npos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot8.pdf", width = 6.5, height = 3)
par(mar = c(1.5,4.5, 2, 9.75))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = sel$Npos, names.arg = NULL, col = sel$metab.color, axes = F, add = T)
barplot(height = -sel$Nneg, names.arg = NULL, col = sel$metab.color, axes = F, add = T)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 18, labels = sel$COG_category, xpd = NA, srt = 0, adj = .5, cex = .7)
ticks <- axis(side = 2, line = .5, labels = F, tck = -.02)
axis(side = 2, line = 0, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes per COG category (%)", side = 2, line = 3.75, cex = .7)
mtext(text = "Positive Selection", side = 2, line = 2.5, at = y.lim[2] / 2, cex = .7)
mtext(text = "Purifying Selection", side = 2, line = 2.5, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .4, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = .25, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 1.25, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 10, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:length(bar.spots)), labels = paste(sel$COG_category,"=",sel$label), 
     xpd = NA, adj = 0, cex = .6, col = sel$metab.color)
dev.off()

# ---- breakdown by season ----

pos <- genes[log.MK < -1.5 & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description, Broad_category, season)]
pos[is.na(COG_category), COG_category := "?"]

neg <- genes[log.MK > 1.5 & is.finite(log.MK), .(.N), by = .(COG_category, COG_Description, Broad_category, season)]
neg[is.na(COG_category), COG_category := "?"]

# OK, what to do about the duplicate (multi-digit) COG categories?
# not sure how to add the multiple-category ones in.. remove for now *************
pos <- pos[nchar(COG_category) == 1]
neg <- neg[nchar(COG_category) == 1]

# And not sure whether to include the unknown functions *******
# ALSO, why are all the NA's not in the S category?
pos <- pos[COG_category != "?" & COG_category != "S"]
neg <- neg[COG_category != "?" & COG_category != "S"]

sel <- merge(x = pos, y = neg, by = c("COG_category","COG_Description","Broad_category","season"), all = T, suffixes = c("pos","neg"))

# ---- plot all COG categories seasonally ----

col.key <- data.table("Broad_category" = c(NA,"Poorly Characterized","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                      "color" = c("grey70","grey70","skyblue3","seagreen","rosybrown3"))
sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)

index <- order(sel$Npos, decreasing = T)
sel <- sel[index]

pos <- dcast(data = sel, formula = COG_category + color + COG_Description ~ season, value.var = "Npos")
key <- pos[ ,1:3]
pos <- as.matrix(pos[ ,-(2:3)], rownames = T)

neg <- dcast(data = sel, formula = COG_category ~ season, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)

all.equal(row.names(pos), row.names(neg))
all.equal(colnames(pos), colnames(neg))
all.equal(row.names(pos), key$COG_category)

bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, axes = F, plot = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(neg, na.rm = T), max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot9.pdf", width = 6.5, height = 3)
par(mar = c(1.25,3.75, 2.75, 9))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, space = c(0,4), ann = F, axisnames = F, col = key$color, add = T)
barplot(height = -neg, beside = T, names.arg = NULL, axes = F, space = c(0,4), ann = F, axisnames = F, col = key$color, add = T)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .3)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes per COG category (%)", side = 2, line = 3, cex = .7, at = 0)
mtext(text = "Positive Selection", side = 2, line = 2, at = 0, cex = .7, adj = 0)
mtext(text = "Purifying Selection", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1.25, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 25, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$COG_Description), 
     xpd = NA, adj = 0, cex = .3)
mtext(text = col.key[-(1:2),Broad_category], col = col.key[-(1:2),color],
      side = 3,line = c(1.2,0.5,-.1), cex = .7, at = x.lim[2] + x.lim[2] / 25, adj = 0)
dev.off()


# ---- plot all metabolism COGs seasonally ----

sel <- sel[Broad_category == "Metabolism"]

metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                                                 "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                                                 "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                           "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                           "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                       "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, catabolism",
                                       "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
sel <- merge(x = sel, y = metab.colors, by = "COG_Description")

pos <- dcast(data = sel, formula = COG_category + metab.color + COG_Description + label ~ season, value.var = "Npos")
key <- pos[ ,1:4]
pos <- as.matrix(pos[ ,-(2:4)], rownames = T)

neg <- dcast(data = sel, formula = COG_category ~ season, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)

all.equal(row.names(pos), row.names(neg))
all.equal(colnames(pos), colnames(neg))
all.equal(row.names(pos), key$COG_category)


bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(neg, na.rm = T), max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot10.pdf", width = 6.5, height = 3)
par(mar = c(1.25,3.75, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
barplot(height = -neg, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .6)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes per COG category (%)", side = 2, line = 3, cex = .7, at = 0)
mtext(text = "Positive Selection", side = 2, line = 2, at = y.lim[2] / 2, cex = .7, adj = .5)
mtext(text = "Purifying Selection", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)
dev.off()


# ---- plot all metabolism COGs seasonally- pos stacked ----

bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F)
x.lim <- c(0, max(bar.spots) + bar.spots[1])
y.lim <- c(0, max(colSums(pos, na.rm = T), na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot11.pdf", width = 6.5, height = 3)
par(mar = c(1.25,1.75, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
text(x = bar.spots, y = 0 - y.lim[2] / 10, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = -1, labels = F, tck = -.02)
axis(side = 2, line = -1.5,  at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Positive Selection (% COG category)", side = 2, line = .75, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 50, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$label), 
     xpd = NA, adj = 0, cex = .7, col = key$metab.color)
dev.off()

# ---- plot all metabolism COGs seasonally- pos beside ----

bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(0, max(pos, na.rm = T))

# ----
pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot12.pdf", width = 6.5, height = 3)
par(mar = c(1.25,2, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
text(x = bar.spots, y = 0 - y.lim[2] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .6)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Positive Selection (% COG category)", side = 2, line = 1.25, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)
dev.off()
