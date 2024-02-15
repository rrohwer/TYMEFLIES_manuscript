# RRR


# ---- libraries and user input ----

library(data.table)
library(lubridate, warn.conflicts = F, quietly = T)

userprefs <- commandArgs(trailingOnly = TRUE)

genes <- userprefs[1] 
anns <- userprefs[2]
abunds <- userprefs[3]
threads <- userprefs[3]
plot.folder <- userprefs[4]

# COMMENT OUT TESTING PATHS
cat("\n\n you left the mac paths defined! \n\n")
genes <- "~/Desktop/test_plots/gene_info_data/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz"
anns <- "~/Desktop/test_plots/gene_info_data/ME2011-09-21_3300043464_group3_bin69.emapper.annotations.gz"
abunds <- "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz"
threads <- "1"
plot.folder <- "~/Desktop/test_plots"

# ---- set up output folders ----

if(!dir.exists(file.path(plot.folder,"Overall_count"))){
  dir.create(file.path(plot.folder,"Overall_count"))
}
if(!dir.exists(file.path(plot.folder,"Overall_count_unknown"))){
  dir.create(file.path(plot.folder,"Overall_count_unknown"))
}
if(!dir.exists(file.path(plot.folder,"Overall_perc"))){
  dir.create(file.path(plot.folder,"Overall_perc"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_count"))){
  dir.create(file.path(plot.folder,"Metabolism_count"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_perc"))){
  dir.create(file.path(plot.folder,"Metabolism_perc"))
}
if(!dir.exists(file.path(plot.folder,"Overall_seasons_count"))){
  dir.create(file.path(plot.folder,"Overall_seasons_count"))
}
if(!dir.exists(file.path(plot.folder,"Overall_seasons_count_unknown"))){
  dir.create(file.path(plot.folder,"Overall_seasons_count_unknown"))
}
if(!dir.exists(file.path(plot.folder,"Overall_seasons_perc"))){
  dir.create(file.path(plot.folder,"Overall_seasons_perc"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_seasons_count"))){
  dir.create(file.path(plot.folder,"Metabolism_seasons_count"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_seasons_perc"))){
  dir.create(file.path(plot.folder,"Metabolism_seasons_perc"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_seasons_pos_count"))){
  dir.create(file.path(plot.folder,"Metabolism_seasons_pos_count"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_seasons_pos_perc"))){
  dir.create(file.path(plot.folder,"Metabolism_seasons_pos_perc"))
}
if(!dir.exists(file.path(plot.folder,"Overall_years_count"))){
  dir.create(file.path(plot.folder,"Overall_years_count"))
}
if(!dir.exists(file.path(plot.folder,"Overall_years_count_unknown"))){
  dir.create(file.path(plot.folder,"Overall_years_count_unknown"))
}
if(!dir.exists(file.path(plot.folder,"Overall_years_perc"))){
  dir.create(file.path(plot.folder,"Overall_years_perc"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_years_count"))){
  dir.create(file.path(plot.folder,"Metabolism_years_count"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_years_perc"))){
  dir.create(file.path(plot.folder,"Metabolism_years_perc"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_years_pos_count"))){
  dir.create(file.path(plot.folder,"Metabolism_years_pos_count"))
}
if(!dir.exists(file.path(plot.folder,"Metabolism_years_pos_perc"))){
  dir.create(file.path(plot.folder,"Metabolism_years_pos_perc"))
}

if(!dir.exists(file.path(plot.folder,"Selection_Summaries"))){
  dir.create(file.path(plot.folder,"Selection_Summaries"))
}


# ---- format data ----

threads <- as.numeric(threads)
genes <- fread(file = genes, nThread = threads)
anns <- fread(file = anns, sep = "\t", quote = "", fill = T, skip = 4, nThread = threads)
anns <- anns[-((nrow(anns) - 2):nrow(anns))] # stupid stats listed at end of file

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
genes[ ,invasion := factor(invasion, levels = c("none","spiny","zebra"))]
genes[ ,sample := sub(pattern = "\\.IS_gene_info\\.tsv", replacement = "", x = sample)]
genes[ ,sample := sub(pattern = ".gz", replacement = "", x = sample)]

abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz", sep = "\t")
colnames(abund)[3] <- "abund"
abund <- as.data.table(abund)
abund <- abund[ ,.(Sample, Genome, abund)]

genes <- merge(x = genes, y = abund, by.x = c("genome","sample"), by.y = c("Genome","Sample"), all = F)

# ---- summarize by COG ----

# group all unclear COG categories together
genes[is.na(COG_category) | COG_category == "" | nchar(COG_category) > 1 | COG_category == "S", `:=`(COG_category = "?", COG_Description = "Unknown Function", Broad_category = "Unknown")] 
anns[COG_category == "-" | COG_category == "S" | nchar(COG_category) > 1, COG_category := "?"]

# tally up selected genes per COG category on each date
pos <- genes[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman), .(.N), by = .(COG_category, COG_Description,Broad_category, sample, abund, date, season, year)]
neg <- genes[MK.p.val >= .05 & mcdonald.kreitman > 1 & is.finite(mcdonald.kreitman), .(.N), by = .(COG_category, COG_Description, Broad_category, sample, abund, date, season, year)]
sel.summary <- merge(x = pos, y = neg, by = c("COG_category","COG_Description","Broad_category", "sample","abund","date", "season", "year"), all = T, suffixes = c("pos","neg"))

# want to say "X percent of this genome's COG category shows evidence of selection in this time period"
tot.per.COG <- anns[ , .(total.in.COG = .N), by = COG_category] # have to use anns file b/c not all genes present in each day of the gene_info files
sel.summary <- merge(x = sel.summary, y = tot.per.COG, by = "COG_category", all.x = T, all.y = F)
sel.summary[ , `:=`(Npos.perc = Npos / total.in.COG * 100, Nneg.perc = Nneg / total.in.COG * 100)]

# further group dates by seasons and years

sel.seasonal <- sel.summary[ , .(Npos = sum(Npos, na.rm = T), Nneg = sum(Nneg, na.rm = T), Npos.perc = sum(Npos.perc, na.rm = T), Nneg.perc = sum(Nneg.perc, na.rm = T)), 
                             by = .(COG_category, COG_Description, Broad_category, season)]

sel.annual <- sel.summary[ , .(Npos = sum(Npos, na.rm = T), Nneg = sum(Nneg, na.rm = T), Npos.perc = sum(Npos.perc, na.rm = T), Nneg.perc = sum(Nneg.perc, na.rm = T)), 
                           by = .(COG_category, COG_Description, Broad_category, year)]

sel.overall <- sel.summary[ , .(Npos = sum(Npos, na.rm = T), Nneg = sum(Nneg, na.rm = T), Npos.perc = sum(Npos.perc, na.rm = T), Nneg.perc = sum(Nneg.perc, na.rm = T)), 
                            by = .(COG_category, COG_Description, Broad_category)]


# ---- plotting functions ----

my.genome <- genes[1, genome]
tax.label <- paste(genes[1, .(phylum, class, order, family, genus, species)], collapse = " ")
plot.label <- paste(tax.label, paste0(my.genome,".pdf"))

plot.abund.overlay <- function(my.sel){
  sel <- copy(my.sel)
  col.key <- data.table(season = factor(x = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")),
                        color = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))
  sel <- merge(x = sel, y = col.key, by = "season")
  
  plot(x = yday(sel$date), y = sel$abund, type = "n", ann = F, axes = F)
  box()
  points(x = yday(sel$date), y = sel$abund, col = adjustcolor(sel$color, alpha.f = .2), pch = 19, cex = 1)
}

plot.overall <- function(my.sel, remove.unknown, use.percent.cog){
  sel <- copy(my.sel)
  if (use.percent.cog){
    sel[ ,`:=`(Npos = Npos.perc, Nneg = Nneg.perc)]
    y.ax.label <- "Genes Under Selection (percent of COG category)"
  }else{
    y.ax.label <- "Genes Under Selection (count)"
  }
  col.key <- data.table("Broad_category" = c("Unknown","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                        "color" = c("grey70","skyblue3","seagreen","rosybrown3"))
  sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)
  index <- order(sel$Npos + sel$Nneg, decreasing = T)
  sel <- sel[index]
  if (remove.unknown){
    sel <- sel[COG_category != "?"]
  }
  bar.spots <- barplot(height = sel$Npos, names.arg = NULL, col = sel$color, axes = F, plot = F)
  x.lim <- c(min(bar.spots), max(bar.spots))
  y.lim <- c(-max(sel$Nneg, na.rm = T), max(sel$Npos, na.rm = T))
  
  par(mar = c(1.5,3.5, 3, 13))
  plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
  barplot(height = sel$Npos, names.arg = NULL, col = sel$color, axes = F, add = T)
  barplot(height = -sel$Nneg, names.arg = NULL, col = sel$color, axes = F, add = T)
  text(x = bar.spots, y = y.lim[1] + y.lim[1] / 15, labels = sel$COG_category, xpd = NA, srt = 0, adj = .5, cex = .7)
  ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
  axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
  mtext(text = y.ax.label, side = 2, line = 2.75, cex = .7)
  mtext(text = "Positive", side = 2, line = 1.75, at = y.lim[2] / 2, cex = .7)
  mtext(text = "Purifying", side = 2, line = 1.75, at = y.lim[1] / 2, cex = .7)
  mtext(text = "COG category", side = 1, line = .4, cex = .7)
  mtext(text = my.genome, side = 3, line = 1.25, adj = 0, cex = .7)
  mtext(text = tax.label, side = 3, line = 2.15, adj = 0, cex = .7)
  text(x = x.lim[2] + x.lim[2] / 20, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:length(bar.spots)), labels = paste(sel$COG_category,"=",sel$COG_Description), 
       xpd = NA, adj = 0, cex = .5)
  mtext(text = col.key[-(1),Broad_category], col = col.key[-(1),color],
        side = 3,line = c(1.2,0.5,-.1), cex = .7, at = x.lim[2] + x.lim[2] / 20, adj = 0)
}

plot.overall.metab <- function(my.sel, use.percent.cog){
  sel <- copy(my.sel)
  sel <- sel[Broad_category == "Metabolism"]
  if (use.percent.cog){
    sel[ ,`:=`(Npos = Npos.perc, Nneg = Nneg.perc)]
    y.ax.label <- "Genes Under Selection (percent of COG category)"
  }else{
    y.ax.label <- "Genes Under Selection (count)"
  }
  metab.colors <- data.table("COG_Description" = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Inorganic ion transport and metabolism",
                                                   "Energy production and conversion", "Lipid transport and metabolism", "Secondary metabolites biosynthesis, transport, and catabolism",
                                                   "Coenzyme transport and metabolism", "Nucleotide transport and metabolism"),
                             "metab.color" = c("aquamarine3", "chocolate2","lightblue3","darkgoldenrod2", "indianred2", "darkmagenta","pink2","royalblue3"),
                             "label" = c("Amino acid\ntransport and metabolism", "Carbohydrate\ntransport and metabolism", "Inorganic ion\ntransport and metabolism",
                                         "Energy\nproduction and conversion", "Lipid\ntransport and metabolism", "Secondary metabolites\nbiosynthesis, transport, and catabolism",
                                         "Coenzyme\ntransport and metabolism", "Nucleotide\ntransport and metabolism"))
  sel <- merge(x = sel, y = metab.colors, by = "COG_Description")
  index <- order(sel$Npos + sel$Nneg, decreasing = T)
  sel <- sel[index]
  bar.spots <- barplot(height = sel$Npos, names.arg = NULL, col = sel$color, axes = F, plot = F)
  x.lim <- c(min(bar.spots), max(bar.spots))
  y.lim <- c(-max(sel$Nneg, na.rm = T), max(sel$Npos, na.rm = T))
  
  par(mar = c(1.5,4.5, 2, 9.75))
  plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
  barplot(height = sel$Npos, names.arg = NULL, col = sel$metab.color, axes = F, add = T)
  barplot(height = -sel$Nneg, names.arg = NULL, col = sel$metab.color, axes = F, add = T)
  text(x = bar.spots, y = y.lim[1] + y.lim[1] / 18, labels = sel$COG_category, xpd = NA, srt = 0, adj = .5, cex = .7)
  ticks <- axis(side = 2, line = .5, labels = F, tck = -.02)
  axis(side = 2, line = 0, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
  mtext(text = y.ax.label, side = 2, line = 3.75, cex = .7)
  mtext(text = "Positive", side = 2, line = 2.5, at = y.lim[2] / 2, cex = .7)
  mtext(text = "Purifying", side = 2, line = 2.5, at = y.lim[1] / 2, cex = .7)
  mtext(text = "COG category", side = 1, line = .4, cex = .7)
  mtext(text = unique(genes$genome), side = 3, line = .25, adj = 0, cex = .7)
  mtext(text = tax.label, side = 3, line = 1.25, adj = 0, cex = .7)
  text(x = x.lim[2] + x.lim[2] / 10, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:length(bar.spots)), labels = paste(sel$COG_category,"=",sel$label), 
       xpd = NA, adj = 0, cex = .6, col = sel$metab.color)
}

plot.seasons <- function(my.sel, remove.unknown, use.percent.cog){
  sel <- copy(my.sel)
  if (use.percent.cog){
    sel[ ,`:=`(Npos = Npos.perc, Nneg = Nneg.perc)]
    y.ax.label <- "Genes Under Selection (percent of COG category)"
  }else{
    y.ax.label <- "Genes Under Selection (count)"
  }
  col.key <- data.table("Broad_category" = c("Unknown","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                        "color" = c("grey70","skyblue3","seagreen","rosybrown3"))
  sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)
  if (remove.unknown){
    sel <- sel[COG_category != "?"]
  }
  pos <- dcast(data = sel, formula = COG_category + color + COG_Description ~ season, value.var = "Npos")
  key <- pos[ ,1:3]
  pos <- as.matrix(pos[ ,-(2:3)], rownames = T)
  neg <- dcast(data = sel, formula = COG_category ~ season, value.var = "Nneg")
  neg <- as.matrix(neg, rownames = T)
  # all.equal(row.names(pos), row.names(neg))
  # all.equal(colnames(pos), colnames(neg))
  # all.equal(row.names(pos), key$COG_category)
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
  mtext(text = y.ax.label, side = 2, line = 3, cex = .7)
  mtext(text = "Positive", side = 2, line = 2, at = 0, cex = .7, adj = 0)
  mtext(text = "Purifying", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
  mtext(text = "COG category", side = 1, line = .25, cex = .7)
  mtext(text = unique(genes$genome), side = 3, line = 1.25, adj = 0, cex = .7)
  mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
  text(x = x.lim[2] + x.lim[2] / 25, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$COG_Description), 
       xpd = NA, adj = 0, cex = .3)
  mtext(text = col.key[-(1),Broad_category], col = col.key[-(1),color],
        side = 3,line = c(1.2,0.5,-.1), cex = .7, at = x.lim[2] + x.lim[2] / 25, adj = 0)
}

plot.seasons.metab <- function(my.sel, use.percent.cog){
  sel <- copy(my.sel)
  if (use.percent.cog){
    sel[ ,`:=`(Npos = Npos.perc, Nneg = Nneg.perc)]
    y.ax.label <- "Genes Under Selection (percent of COG category)"
  }else{
    y.ax.label <- "Genes Under Selection (count)"
  }
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
  mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
  mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
  text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
       xpd = NA, adj = 0, cex = .6, col = key$metab.color)
}

plot.seasons.metab.pos <- function(my.sel, use.percent.cog){
  sel <- copy(my.sel)
  if (use.percent.cog){
    sel[ ,`:=`(Npos = Npos.perc, Nneg = Nneg.perc)]
    y.ax.label <- "Genes Under Positive Selection (% COG category)"
  }else{
    y.ax.label <- "Genes Under Positive Selection (count)"
  }
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
  bar.spots <- barplot(height = pos, beside = T, names.arg = NULL, col = key$metab.color, axes = F, plot = F, axisnames = F, space = c(0,4))
  x.lim <- c(min(bar.spots), max(bar.spots))
  y.lim <- c(0, max(pos, na.rm = T))
  
  par(mar = c(1.25,2, 2.75, 8.25))
  plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
  barplot(height = pos, beside = T, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T, space = c(0,4))
  text(x = bar.spots, y = 0 - y.lim[2] / 20, labels = key$COG_category, xpd = NA, srt = 0, adj = .5, cex = .6)
  text(x = bar.spots[floor(nrow(bar.spots)/2), ], y = y.lim[2] + y.lim[2] / 15, labels = colnames(pos), cex = .7, xpd = NA)
  ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
  axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
  mtext(text = y.ax.label, side = 2, line = 1.25, cex = .7)
  mtext(text = "COG category", side = 1, line = .25, cex = .7)
  mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
  mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
  text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label),
       xpd = NA, adj = 0, cex = .6, col = key$metab.color)
}

plot.years <- function(my.sel, remove.unknown, use.percent.cog){
  sel <- copy(my.sel)
  if (use.percent.cog){
    sel[ ,`:=`(Npos = Npos.perc, Nneg = Nneg.perc)]
    y.ax.label <- "Genes Under Selection (percent of COG category)"
  }else{
    y.ax.label <- "Genes Under Selection (count)"
  }
  col.key <- data.table("Broad_category" = c("Unknown","Information Storage and Processing","Cellular Processes and Signaling","Metabolism"),
                        "color" = c("grey70","skyblue3","seagreen","rosybrown3"))
  sel <- merge(x = sel, y = col.key, by = "Broad_category", all.x = T, all.y = F)
  if (remove.unknown){
    sel <- sel[COG_category != "?"]
  }
  pos <- dcast(data = sel, formula = COG_category + color + COG_Description ~ year, value.var = "Npos")
  pos <- pos[order(color)]
  key <- pos[ ,1:3]
  pos <- as.matrix(pos[ ,-(2:3)], rownames = T)
  pos[is.na(pos)]<- 0
  neg <- dcast(data = sel, formula = COG_category + color + COG_Description ~ year, value.var = "Nneg")
  neg <- neg[order(color)]
  neg <- as.matrix(neg[ ,-(2:3)], rownames = T)
  neg[is.na(neg)] <- 0
  # all.equal(row.names(pos), row.names(neg))
  # all.equal(key$COG_category, row.names(pos))
  # all.equal(colnames(pos), colnames(neg))
  bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, axes = F, plot = F)
  x.lim <- c(min(bar.spots), max(bar.spots))
  y.lim <- c(-max(colSums(neg, na.rm = T)), max(colSums(pos, na.rm = T)))
  
  par(mar = c(1.25,3.75, 2.75, 1))
  plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
  barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$color, add = T)
  barplot(height = -neg, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$color, add = T)
  abline(h = 0)
  text(x = bar.spots, y = y.lim[1] + y.lim[1] / 10, labels = colnames(pos), cex = .7, xpd = NA, srt = 90)
  ticks <- axis(side = 2, at = c(-1200,-800,-400,0,200), line = 0, labels = F, tck = -.02)
  axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
  mtext(text = y.ax.label, side = 2, line = 3, cex = .7)
  mtext(text = "Positive", side = 2, line = 2, at = 0, cex = .7, adj = 0)
  mtext(text = "Purifying", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
  mtext(text = unique(genes$genome), side = 3, line = 1.25, adj = 0, cex = .7)
  mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
  mtext(text = col.key[-(1),Broad_category], col = col.key[-(1),color],
        side = 3,line = c(1.25,0.55,-.05), cex = .7, at = x.lim[2] - x.lim[2] / 4, adj = 0)
}

plot.years.metab <- function(my.sel, use.percent.cog){
  sel <- copy(my.sel)
  if (use.percent.cog){
    sel[ ,`:=`(Npos = Npos.perc, Nneg = Nneg.perc)]
    y.ax.label <- "Genes Under Selection (percent of COG category)"
  }else{
    y.ax.label <- "Genes Under Selection (count)"
  }
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
  pos[is.na(pos)] <- 0
  neg <- dcast(data = sel, formula = COG_category ~ year, value.var = "Nneg")
  neg <- as.matrix(neg, rownames = T)
  neg[is.na(neg)] <- 0
  # all.equal(row.names(pos), row.names(neg))
  # all.equal(colnames(pos), colnames(neg))
  # all.equal(row.names(pos), key$COG_category)
  bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, col = key$metab.color, axes = F, plot = F, axisnames = F)
  x.lim <- c(min(bar.spots), max(bar.spots))
  y.lim <- c(-max(colSums(neg, na.rm = T)), max(colSums(pos, na.rm = T)))
  
  par(mar = c(1.25,3.75, 2.75, 8.25))
  plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
  barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
  barplot(height = -neg, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
  abline(h = 0)
  text(x = bar.spots, y = y.lim[1] + y.lim[1] / 10, labels = colnames(pos), cex = .7, xpd = NA, srt = 90)
  ticks <- axis(side = 2, at = c(-1200,-800,-400,0,200), line = 0, labels = F, tck = -.02)
  axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
  mtext(text = y.ax.label, side = 2, line = 3, cex = .7)
  mtext(text = "Positive", side = 2, line = 2, at = y.lim[2] / 2, cex = .7, adj = .5)
  mtext(text = "Purifying", side = 2, line = 2, at = y.lim[1] / 2, cex = .7)
  mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
  mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
  text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label), 
       xpd = NA, adj = 0, cex = .6, col = key$metab.color)
}

plot.years.metab.pos <- function(my.sel, use.percent.cog){
  sel <- copy(my.sel)
  if (use.percent.cog){
    sel[ ,`:=`(Npos = Npos.perc, Nneg = Nneg.perc)]
    y.ax.label <- "Genes Under Positive Selection (% of COG category)"
  }else{
    y.ax.label <- "Genes Under Positive Selection (count)"
  }
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
  pos[is.na(pos)] <- 0
  bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, col = key$metab.color, axes = F, plot = F, axisnames = F)
  x.lim <- c(min(bar.spots), max(bar.spots))
  y.lim <- c(0, max(colSums(pos, na.rm = T)))
  
  par(mar = c(1.25,2, 2.75, 8.25))
  plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
  barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
  text(x = bar.spots, y = 0 - y.lim[2] / 20, labels = colnames(pos), cex = .7, xpd = NA, srt = 90)
  ticks <- axis(side = 2, line = 0, labels = F, tck = -.02)
  axis(side = 2, line = -.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
  mtext(text = y.ax.label, side = 2, line = 1.25, cex = .7)
  mtext(text = "COG category", side = 1, line = .25, cex = .7)
  mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
  mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
  text(x = x.lim[2] + x.lim[2] / 15, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$COG_category,"=",key$label),
       xpd = NA, adj = 0, cex = .6, col = key$metab.color)
}


# ---- export plots ----

pdf(file = file.path(plot.folder,"Overall_count",plot.label), width = 6.5, height = 3)
par(fig = c(0,1,.2,1))
plot.overall(my.sel = sel.overall, remove.unknown = T, use.percent.cog = F)
par(fig = c(0,1,0,.2), mar = c(1,1,1,1), new = T)
plot.abund.overlay(my.sel = sel.summary)
dev.off()

pdf(file = file.path(plot.folder,"Overall_count_unknown",plot.label), width = 6.5, height = 3)
plot.overall(my.sel = sel.overall, remove.unknown = F, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Overall_perc",plot.label), width = 6.5, height = 3)
plot.overall(my.sel = sel.overall, remove.unknown = F, use.percent.cog = T)
dev.off()

pdf(file = file.path(plot.folder,"Metabolism_count", plot.label), width = 6.5, height = 3)
plot.overall.metab(my.sel = sel.overall, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Metabolism_perc",plot.label), width = 6.5, height = 3)
plot.overall.metab(my.sel = sel.overall, use.percent.cog = T)
dev.off()

pdf(file = file.path(plot.folder,"Overall_seasons_count", plot.label), width = 6.5, height = 3)
plot.seasons(my.sel = sel.seasonal, remove.unknown = T, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Overall_seasons_count_unknown", plot.label), width = 6.5, height = 3)
plot.seasons(my.sel = sel.seasonal, remove.unknown = F, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Overall_seasons_perc",plot.label), width = 6.5, height = 3)
plot.seasons(my.sel = sel.seasonal, remove.unknown = F, use.percent.cog = T)
dev.off()

pdf(file = file.path(plot.folder,"Metabolism_seasons_count",plot.label), width = 6.5, height = 3)
plot.seasons.metab(my.sel = sel.seasonal, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Metabolism_seasons_perc",plot.label), width = 6.5, height = 3)
plot.seasons.metab(my.sel = sel.seasonal, use.percent.cog = T)
dev.off()

pdf(file = file.path(plot.folder,"Metabolism_seasons_pos_count",plot.label), width = 6.5, height = 3)
plot.seasons.metab.pos(my.sel = sel.seasonal, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Metabolism_seasons_pos_perc",plot.label), width = 6.5, height = 3)
plot.seasons.metab.pos(my.sel = sel.seasonal, use.percent.cog = T)
dev.off()

pdf(file = file.path(plot.folder,"Overall_years_count", plot.label), width = 6.5, height = 3)
plot.years(my.sel = sel.annual, remove.unknown = T, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Overall_years_count_unknown", plot.label), width = 6.5, height = 3)
plot.years(my.sel = sel.annual, remove.unknown = F, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Overall_years_perc",plot.label), width = 6.5, height = 3)
plot.years(my.sel = sel.annual, remove.unknown = F, use.percent.cog = T)
dev.off()

pdf(file = file.path(plot.folder,"Metabolism_years_count",plot.label), width = 6.5, height = 3)
plot.years.metab(my.sel = sel.annual, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Metabolism_years_perc",plot.label), width = 6.5, height = 3)
plot.years.metab(my.sel = sel.annual, use.percent.cog = T)
dev.off()

pdf(file = file.path(plot.folder,"Metabolism_years_pos_count",plot.label), width = 6.5, height = 3)
plot.years.metab.pos(my.sel = sel.annual, use.percent.cog = F)
dev.off()
pdf(file = file.path(plot.folder,"Metabolism_years_pos_perc",plot.label), width = 6.5, height = 3)
plot.years.metab.pos(my.sel = sel.annual, use.percent.cog = T)
dev.off()

# # ---- plot all metabolism COGs seasonally- pos stacked ----
#
# bar.spots <- barplot(height = pos, beside = F, names.arg = NULL, col = key$metab.color, axes = F, plot = T, axisnames = F)
# x.lim <- c(0, max(bar.spots) + bar.spots[1])
# y.lim <- c(0, max(colSums(pos, na.rm = T), na.rm = T))
#
# # ----
# pdf(file = "data/2023-04-19_test_merging_gene_info_with_annotations/test_plot5.pdf", width = 6.5, height = 3)
# par(mar = c(1.25,1.75, 2.75, 8.25))
# plot(x = 1, y = 1, xlim = x.lim, ylim = y.lim, type = "n", axes = F, ann = F)
# barplot(height = pos, beside = F, names.arg = NULL, axes = F, ann = F, axisnames = F, col = key$metab.color, add = T)
# text(x = bar.spots, y = 0 - y.lim[2] / 10, labels = colnames(pos), cex = .7, xpd = NA)
# ticks <- axis(side = 2, line = -1, labels = F, tck = -.02)
# axis(side = 2, line = -1.5, at = ticks, labels = abs(ticks), lwd = 0, las = 2, cex.axis = .7)
# mtext(text = "Genes Under Positive Selection (count)", side = 2, line = .75, cex = .7)
# mtext(text = unique(genes$genome), side = 3, line = 1, adj = 0, cex = .7)
# mtext(text = tax.label, side = 3, line = 2, adj = 0, cex = .7)
# text(x = x.lim[2] + x.lim[2] / 50, y = seq(from = y.lim[2], to = y.lim[1], along.with = 1:nrow(key)), labels = paste(key$label),
#      xpd = NA, adj = 0, cex = .7, col = key$metab.color)
# dev.off()

# ---- export data ----

sel.summary[ ,date := as.character(date)]

fwrite(x = sel.summary, file = file.path(plot.folder,"Selection_Summaries", paste0(my.genome,"_selected_COGs.tsv.gz")), sep = "\t", row.names = F, nThread = threads, compress = "gzip")

cat("Finished:", my.genome, "\n")

# ----