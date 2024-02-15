# RRR
# are there phylogenetic or seasonal trends in which COGs get selected more?

library(data.table)

# ---- import ----

# get tax
tax <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
tax <- as.data.table(tax)
tax <- tax[winner == TRUE, .(bin.full.name, domain, phylum, class, order, family, genus, species)]

# get COGs
folder.path <- "data/2023-05-31_COG_summary_data/Selection_Summaries/"
my.files <- list.files(folder.path)

all.genomes.list <- list()
for (f in my.files){
  one.g <- fread(file = file.path(folder.path,f))
  if (nrow(one.g > 0)){ # some had no selection so are empty tables
    one.g$genome <- sub(pattern = "_selected_COGs\\.tsv\\.gz", replacement = "", x = f)
    all.genomes.list <- c(all.genomes.list, list(one.g))
  }
}

all.genomes.table <- rbindlist(l = all.genomes.list)
# rm(all.genomes.list)
all.genomes.table <- merge(x = all.genomes.table, y = tax, by.x = "genome", by.y = "bin.full.name", all = T)

# ---- overall seasonal trends? ----

seas <- all.genomes.table[ , .(Npos = sum(Npos, na.rm = T), Nneg = sum(Nneg, na.rm = T), Npos.perc = sum(Npos.perc, na.rm = T), Nneg.perc = sum(Nneg.perc, na.rm = T)), by = .(season, COG_category, COG_Description, Broad_category)]

ann <- all.genomes.table[ , .(Npos = sum(Npos, na.rm = T), Nneg = sum(Nneg, na.rm = T), Npos.perc = sum(Npos.perc, na.rm = T), Nneg.perc = sum(Nneg.perc, na.rm = T)), by = .(year, COG_category, COG_Description, Broad_category)]

phyla <- all.genomes.table[ , .(Npos = sum(Npos, na.rm = T), Nneg = sum(Nneg, na.rm = T), Npos.perc = sum(Npos.perc, na.rm = T), Nneg.perc = sum(Nneg.perc, na.rm = T)), by = .(phylum, COG_category, COG_Description, Broad_category)]

# how do I normalize the phylum numbers???

# ---- quick looks - phylum trends? -----

sel <- copy(phyla)
col.key <- data.table("Broad_category" = c("Unknown","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                      "color" = c("grey70","skyblue3","seagreen","rosybrown3"))
sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)
pos <- dcast(data = sel, formula = COG_category + color + COG_Description ~ phylum, value.var = "Npos")
pos <- pos[order(color)]
key <- pos[ ,1:3]
pos <- as.matrix(pos[ ,-(2:3)], rownames = T)
pos[is.na(pos)]<- 0
neg <- dcast(data = sel, formula = COG_category + color + COG_Description ~ phylum, value.var = "Nneg")
neg <- neg[order(color)]
neg <- as.matrix(neg[ ,-(2:3)], rownames = T)
neg[is.na(neg)] <- 0
bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, axes = F, plot = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(colSums(neg, na.rm = T)), max(colSums(pos, na.rm = T)))

par(mar = c(5,3.75, 2.75, 1))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$color, add = T)
barplot(height = -neg, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$color, add = T)
abline(h = 0)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 10, labels = colnames(pos), cex = .7, xpd = NA, srt = 90)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text =  "Genes Under Selection (count)", side = 2, line = 3, cex = .7)
mtext(text = "Positive", side = 2, line = 2, at = 0, cex = .7, adj = 0)
mtext(text = "Purifying", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = col.key[-(1),Broad_category], col = col.key[-(1),color],
      side = 3,line = c(1.25,0.55,-.05), cex = .7, at = x.lim[2] - x.lim[2] / 4, adj = 0)

# seems like the abundant phyla have more selection, probably because they're abundant enough to see it?
# not sure how to normalize to see if some phyla have MORE selection

sel <- copy(phyla)
sel <- sel[Broad_category == "Metabolism"]
metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                                                 "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                                                 "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                           "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                           "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                       "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, catabolism",
                                       "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
sel <- merge(x = sel, y = metab.colors, by = "COG_Description")
pos <- dcast(data = sel, formula = COG_category + metab.color + COG_Description + label ~ phylum, value.var = "Npos")
key <- pos[ ,1:4]
pos <- as.matrix(pos[ ,-(2:4)], rownames = T)
pos[is.na(pos)] <- 0
neg <- dcast(data = sel, formula = COG_category ~ phylum, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)
neg[is.na(neg)] <- 0
bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, col = key$metab.color, axes = F, plot = F, axisnames = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(colSums(neg, na.rm = T)), max(colSums(pos, na.rm = T)))

par(mar = c(5,3.75, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
barplot(height = -neg, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
abline(h = 0)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 10, labels = colnames(pos), cex = .7, xpd = NA, srt = 90)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Selection", side = 2, line = 3, cex = .7)
mtext(text = "Positive", side = 2, line = 2, at = y.lim[2] / 2, cex = .7, adj = .5)
mtext(text = "Purifying", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)




sel <- copy(phyla)
sel <- sel[Broad_category == "Metabolism"]
metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                                                 "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                                                 "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                           "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                           "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                       "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, catabolism",
                                       "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
sel <- merge(x = sel, y = metab.colors, by = "COG_Description")
pos <- dcast(data = sel, formula = COG_category + metab.color + COG_Description + label ~ phylum, value.var = "Npos")
key <- pos[ ,1:4]
pos <- as.matrix(pos[ ,-(2:4)], rownames = T)
pos[is.na(pos)] <- 0
bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, col = key$metab.color, axes = F, plot = F, axisnames = F)
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(0, max(colSums(pos, na.rm = T)))

par(mar = c(5,2, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
text(x = bar.spots, y = 0 - y.lim[2] / 10, labels = colnames(pos), cex = .7, xpd = NA, srt = 90)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Positive Selection (count)", side = 2, line = 1.25, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label),
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)

# ---- quick looks- season trends? ----

sel <- copy(seas)
col.key <- data.table("Broad_category" = c("Unknown","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                      "color" = c("grey70","skyblue3","seagreen","rosybrown3"))
sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)
pos <- dcast(data = sel, formula = COG_category + color + COG_Description ~ season, value.var = "Npos")
key <- pos[ ,1:3]
pos <- as.matrix(pos[ ,-(2:3)], rownames = T)
neg <- dcast(data = sel, formula = COG_category ~ season, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)
bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, axes = F, plot = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(neg, na.rm = T), max(pos, na.rm = T))

par(mar = c(1.25,3.75, 2.75, 9))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, space = c(0,4), ann = F, axisnames = F, col = key$color, add = T)
barplot(height = -neg, beside = T, names.arg = NULL, axes = F, space = c(0,4), ann = F, axisnames = F, col = key$color, add = T)
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .3)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = "Genes Under Selection (count)", side = 2, line = 3, cex = .7)
mtext(text = "Positive", side = 2, line = 2, at = 0, cex = .7, adj = 0)
mtext(text = "Purifying", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
text(x = x.lim[2] + x.lim[2] / 25, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$COG_Description), 
     xpd = NA, adj = 0, cex = .3)
mtext(text = col.key[-(1),Broad_category], col = col.key[-(1),color],
      side = 3,line = c(1.2,0.5,-.1), cex = .7, at = x.lim[2] + x.lim[2] / 25, adj = 0)


sel <- copy(seas)
sel <- sel[Broad_category == "Metabolism"]
metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                                                 "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                                                 "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                           "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                           "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                       "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, catabolism",
                                       "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
sel <- merge(x = sel, y = metab.colors, by = "COG_Description")
sel[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
sel <- sel[order(season)]
pos <- dcast(data = sel, formula = COG_category + metab.color + COG_Description + label ~ season, value.var = "Npos")
key <- pos[ ,1:4]
pos <- as.matrix(pos[ ,-(2:4)], rownames = T)
neg <- dcast(data = sel, formula = COG_category ~ season, value.var = "Nneg")
neg <- as.matrix(neg, rownames = T)
# all.equal(row.names(pos), row.names(neg))
# all.equal(colnames(pos), colnames(neg))
# all.equal(row.names(pos), key$COG_category)
bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, col = key$metab.color, axes = F, plot = F, axisnames = F, space = c(0,4))
x.lim <- c(min(bar.spots), max(bar.spots))
y.lim <- c(-max(neg, na.rm = T), max(pos, na.rm = T))

par(mar = c(1.25,3.75, 2.75, 8.25))
plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
barplot(height = pos, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
barplot(height = -neg, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
text(x = bar.spots, y = y.lim[1] + y.lim[1] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .6)
text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 3, labels = colnames(pos), cex = .7, xpd = NA)
ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
mtext(text = y.ax.label, side = 2, line = 3, cex = .7)
mtext(text = "Positive", side = 2, line = 2, at = y.lim[2] / 2, cex = .7, adj = .5)
mtext(text = "Purifying", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
mtext(text = "COG category", side = 1, line = .25, cex = .7)
mtext(text = "Genes Under Selection (count)", side = 3, line = 2, adj = 0, cex = .7)
text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
     xpd = NA, adj = 0, cex = .6, col = key$metab.color)


# hmmm, not optimistic about my "tracks growth" idea, seems like maybe it tracks diversity

