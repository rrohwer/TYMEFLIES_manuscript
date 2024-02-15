# RRR

# installing ggplot2 and ggtree would be a major pain on the servers
# plus paneling is limited down the road by patchworks alignments
# plus apparently the base heatmap() function is good: https://r-graph-gallery.com/heatmap


# ---- set-up ----

userinput <- commandArgs(trailingOnly = TRUE)
genes.file <- userinput[1]
threads <- as.numeric(userinput[2])
plot.folder <- userinput[3]

library(data.table, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)

# # local path test
# genes.file <- "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.consistently_selected_genes-by_sample.tsv.gz"
# threads = 1
# plot.folder <- "figures/2023-10-11_heatmaps_test_plots"
# 
# genes.file <- "data/2023-09-13_KEGG_annotations/example_files/ME2013-06-27_3300042327_group4_bin82.consistently_selected_genes-by_sample.tsv.gz"

genome <- sub("^.*/","",genes.file)
genome <- sub("\\..*$","",genome)
cat("\nworking on:",genome,"\n")

genes <- fread(file = genes.file, sep = "\t", nThread = threads)

taxon <- paste(unique(genes[ ,.(phylum,class,order,family,genus,species,num.in.cluster)]), collapse = " ")

# ---- define functions ----

make.plot.folder <- function(subfolder, plot.folder){
  if (!dir.exists(file.path(plot.folder,subfolder))){
    cat("making",file.path(plot.folder,subfolder),"\n")
    dir.create(file.path(plot.folder,subfolder))
  }
}

order.genes.by.appearance.pattern <- function(genes.mat){
  genes.dist <- dist(x = genes.mat, method = "euclidean", diag = T)
  genes.clust <- hclust(d = genes.dist, method = "ward.D")
  genes.dend <- as.dendrogram(genes.clust)
  index <- order.dendrogram(genes.dend)
  return(index)
}

order.samples.by.date <- function(genes.mat){
  sample.dates <- colnames(genes.mat) |>
    substr(start = 3, stop = 12) |>
    parse_date_time(orders = "ymd")
  index <- order(sample.dates)
  return(index)
}

order.samples.by.season <- function(genes.mat){
  # order by date first to ensure will be in order within a season
  season <- sub("^.*_", replacement = "", x = colnames(genes.mat)) |>
    factor(levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))
  index <- order(season)
  return(index)
}

get.x.year.ticks <- function(genes.mat){
  sample.year <- colnames(genes.mat) |>
    substr(start = 3, stop = 12) |>
    parse_date_time(orders = "ymd") |>
    year()
  tic.locs <- which(!duplicated(sample.year))
  tic.labs <- sub(pattern = "\\.", replacement = "\n", sample.year[tic.locs])
  tic.labs <- sub("Ice\nOn", "Ice-On", x = tic.labs)
  tics.plus.end <- c(tic.locs, length(sample.year))
  lab.locs <- ((tics.plus.end[-1] - tics.plus.end[-length(tics.plus.end)] ) / 2) + tics.plus.end[-length(tics.plus.end)]
  return(list("tick.locs" = tic.locs, "lab.locs" = lab.locs, "labs" = tic.labs))
}

get.x.season.ticks <- function(genes.mat){
  season <- sub("^.*_", replacement = "", x = colnames(genes.mat)) |>
    factor(levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))
  tic.locs <- which(!duplicated(season))
  tic.labs <- sub(pattern = "\\.", replacement = "\n", season[tic.locs])
  tic.labs <- sub("Ice\nOn", "Ice-On", x = tic.labs)
  tics.plus.end <- c(tic.locs, length(season))
  lab.locs <- ((tics.plus.end[-1] - tics.plus.end[-length(tics.plus.end)] ) / 2) + tics.plus.end[-length(tics.plus.end)]
  return(list("tick.locs" = tic.locs, "lab.locs" = lab.locs, "labs" = tic.labs))
}

make.timeline.plot <- function(genes.mat, col.key, x.axis){
  image.mat <- t(genes.mat)
  image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
  gene.lab <- sub(pattern = genome, replacement = "", x = colnames(image.mat))
  gene.lab <- sub(pattern = "^_scaffold_", replacement = "", x = gene.lab)
  if (length(x.axis$labs) == 20){
    x.axis$labs <- c("2000",NA,"2002",NA,"2004",NA,"2006",NA,"2008",NA,"2010",NA,"2012","2013","2014","2015","2016","2017","2018","2019")
  }
  par(fig = c(0,.905,0,1), mar = c(1.7,3.5,1.6,0))
  image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F,
        col = col.key$color)
  box(which = "plot", lwd = 1)
  axis(side = 1, at = x.axis$tick.locs, labels = F, lwd.ticks = 1, lwd = 0, tck = -.015)
  axis(side = 1, at = x.axis$tick.locs, labels = x.axis$labs, lwd = 0, las = 2, hadj = .5, line = 0, cex.axis = .7)
  axis(side = 2, at = 1:nrow(genes.mat), lwd = 0, labels = gene.lab, las = 2, line = -.75, cex.axis = .4)
  mtext(text = "Consistently Selected Genes", side = 2, line = 2.75, cex = .7)
  mtext(text = taxon, side = 3, adj = 0, cex = .5, line = .2)
  mtext(text = genome, side = 3, adj = 0, cex = .7, line = .8)
  
  par(fig = c(.895,1,.3,.7), mar = c(1,1,1,2), new = T)
  plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
  segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
  axis(side = 4, at = c(.05,.15,.25), xpd = T, lwd = 0, las = 2, line = -.75, cex.axis = .7)
  axis(side = 4, at = c(.05,.15,.25), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
  mtext(text = "Positive\nSelection\np-value", adj = 0, line = 1, cex = .7)
}

make.season.plot <- function(genes.mat, col.key, x.axis){
  image.mat <- t(genes.mat)
  image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this flipping around way it matches the table key 
  gene.lab <- sub(pattern = genome, replacement = "", x = colnames(image.mat))
  gene.lab <- sub(pattern = "^_scaffold_", replacement = "", x = gene.lab)
  par(fig = c(0,.905,0,1), mar = c(3,3.55,1.6,0))
  image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F,
        col = col.key$color)
  box(which = "plot", lwd = 1)
  axis(side = 1, at = x.axis$tick.locs, labels = F, lwd.ticks = 1, lwd = 0, tck = 1)
  axis(side = 1, at = x.axis$tick.locs, labels = F, lwd.ticks = 1, lwd = 0, tck = -.05)
  axis(side = 1, at = x.axis$lab.locs, labels = x.axis$labs, lwd = 0, las = 2, hadj = 1, line = -.7, cex.axis = .7)
  axis(side = 2, at = 1:nrow(genes.mat), lwd = 0, labels = gene.lab, las = 2, line = -.75, cex.axis = .4)
  mtext(text = "Consistently Selected Genes", side = 2, line = 2.75, cex = .7)
  mtext(text = taxon, side = 3, adj = 0, cex = .5, line = .2)
  mtext(text = genome, side = 3, adj = 0, cex = .7, line = .8)
  
  par(fig = c(.895,1,.3,.7), mar = c(1,1,1,2), new = T)
  plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
  segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
  axis(side = 4, at = c(.05,.15,.25), xpd = T, lwd = 0, las = 2, line = -.5, cex.axis = .7)
  axis(side = 4, at = c(.05,.15,.25), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.5, line = 0, labels = F, cex.axis = .7)
  mtext(text = "Positive\nSelection\np-value", adj = 0, line = 1, cex = .7)
}

# ---- format data ----

# set up plot folders, if don't exist
if (!dir.exists(plot.folder)){
  cat("making",plot.folder,"\n")
  dir.create(plot.folder)
}
make.plot.folder(subfolder = "overall-timeline", plot.folder = plot.folder)
make.plot.folder(subfolder = "overall-seasons", plot.folder = plot.folder)
make.plot.folder(subfolder = "2012-timeline", plot.folder = plot.folder)
make.plot.folder(subfolder = "seasons-seasons", plot.folder = plot.folder)

# format table
genes[ ,date := parse_date_time(x = date, orders = "ymd")]

# only call it pos selection if p-value is <= .25, otherwise p-value is 1
genes$positive.selection.pvalue <- as.numeric(1)
genes[mcdonald.kreitman < 1 & MK.p.val <= .25, positive.selection.pvalue := MK.p.val]

# ---- subset to overall ----

genes.group <- genes[consistent.in == "overall"]

# sort genes based on occurrence patterns over time
if(nrow(genes.group) > 0){
  genes.wide <- dcast(data = genes.group, formula = gene + ko + ko.description + module + module.description + path + pathway.description ~ sample + season, value.var = "positive.selection.pvalue")
  genes.mat <- as.matrix(x = genes.wide[ ,-c(2:7)], rownames = 1)
  genes.key <- genes.wide[ ,1:7]
  
  if (nrow(genes.mat) >= 2){
    index <- order.genes.by.appearance.pattern(genes.mat = genes.mat)
    genes.mat <- genes.mat[index, , drop = F]
    genes.key <- genes.key[index]
  }
  
  # go back from p-val = 1 to high p-vals are NA
  genes.mat[genes.mat == 1] <- NA
  
  # make timeline plot ----
  
  # sort dates
  index <- order.samples.by.date(genes.mat = genes.mat)
  genes.mat <- genes.mat[ ,index, drop = F]
  
  # get plot details
  x.axis <- get.x.year.ticks(genes.mat = genes.mat)
  
  n.colors <- 100
  color.fun <- colorRampPalette(colors = c("magenta2", "grey","white"), bias = 1.5)
  col.key <- data.table("color" = color.fun(n = n.colors),
                        "value" = seq(from = min(genes.mat, na.rm = T), to = max(genes.mat, na.rm = T), along.with = 1:n.colors))
  
  # export plot and paired gene key
  if(nrow(genes.mat) > 0){
    pdf(file = file.path(plot.folder,"overall-timeline",paste0(taxon," ",genome,".pdf")), width = 6.5, height = 4)
    make.timeline.plot(genes.mat = genes.mat, col.key = col.key, x.axis = x.axis)
    dev.off()
  }
  
  fwrite(file = file.path(plot.folder,"overall-timeline",paste0(taxon," ",genome,".csv")), x = genes.key, sep = ",", quote = TRUE)
  
  # make seasons plot ----
  
  # sort dates
  index <- order.samples.by.date(genes.mat = genes.mat)
  genes.mat <- genes.mat[ ,index, drop = F]
  
  index <- order.samples.by.season(genes.mat = genes.mat)
  genes.mat <- genes.mat[ ,index, drop = F]
  
  # get plot details
  x.axis <- get.x.season.ticks(genes.mat = genes.mat)
  
  n.colors <- 100
  color.fun <- colorRampPalette(colors = c("magenta2", "grey","white"), bias = 1.5)
  col.key <- data.table("color" = color.fun(n = n.colors),
                        "value" = seq(from = min(genes.mat, na.rm = T), to = max(genes.mat, na.rm = T), along.with = 1:n.colors))
  
  # export plot and paired gene key
  
  if(nrow(genes.mat > 0)){
    pdf(file = file.path(plot.folder,"overall-seasons",paste0(taxon," ",genome,".pdf")), width = 6.5, height = 4)
    make.season.plot(genes.mat = genes.mat, col.key = col.key, x.axis = x.axis)
    dev.off()  
  }
  
  fwrite(file = file.path(plot.folder,"overall-seasons",paste0(taxon," ",genome,".csv")), x = genes.key, sep = ",", quote = TRUE)
  
}


# ---- subset to 2012-groups ----

genes.group <- genes[consistent.in == "pre-2012" | consistent.in == "post-2012"]

if(nrow(genes.group) > 0){
  # sort genes based on occurrence patterns over time (this time need to deal with duplicates from consistent in >1 group)
  genes.wide <- dcast(data = genes.group, formula = gene + consistent.in + ko + ko.description + module + module.description + path + pathway.description ~ sample + season, value.var = "positive.selection.pvalue")
  genes.wide <- genes.wide[order(consistent.in)]
  genes.key <- genes.wide[ ,.(consistent.in = paste(consistent.in, collapse = "; ")), by = .(gene,ko,ko.description,module,module.description,path,pathway.description)]
  genes.mat <- as.matrix(x = genes.wide[ ,-c(2:8)], rownames = 1)
  genes.mat <- genes.mat[!duplicated(row.names(genes.mat)), , drop = F]
  if(!all.equal(row.names(genes.mat), genes.key$gene)){cat("\nTHERE'S A BUG IN LINE 244- key order != matrix order\n")}
  
  # re-order the 2012 groups
  genes.key <- genes.key[ ,consistent.in := factor(x = consistent.in, levels = c("pre-2012","post-2012; pre-2012","post-2012"), ordered = TRUE)]
  index <- order(genes.key$consistent.in)
  genes.key <- genes.key[index, ]
  genes.mat <- genes.mat[index, , drop = F]
  if(!all.equal(row.names(genes.mat), genes.key$gene)){cat("\nTHERE'S A BUG IN LINE 251- key order != matrix order\n")}
  
  # sort genes by appearance, within each group
  my.groups <- unique(genes.key$consistent.in)
  genes.key$old.order <- 1:nrow(genes.key)
  genes.list <- list()
  for (g in my.groups){
    genes.g <- genes.key[consistent.in == g]
    if (nrow(genes.g) >= 2){
      genes.g$new.order <- order.genes.by.appearance.pattern(genes.mat = genes.mat[genes.g$old.order, , drop = F])
    }else{
      genes.g$new.order <- 1
    }
    genes.list <- c(genes.list, list(genes.g))
  }
  genes.list <- rbindlist(l = genes.list)
  
  # re-order matrix and gene key to be ordered by group & gene appearance
  genes.list <- genes.list[order(genes.list$consistent.in, genes.list$new.order)]
  genes.list$new.order <- 1:nrow(genes.list)
  genes.mat <- genes.mat[genes.list$new.order, , drop = F]
  genes.key <- genes.key[genes.list$new.order, ]
  if(!all.equal(row.names(genes.mat), genes.key$gene)){cat("\nTHERE'S A BUG IN LINE 273- key order != matrix order\n")}
  
  # go back from p-val = 1 to high p-vals are NA
  genes.mat[genes.mat == 1] <- NA
  
  # make timeline plot ----
  
  # sort dates
  index <- order.samples.by.date(genes.mat = genes.mat)
  genes.mat <- genes.mat[ ,index, drop = F] 
  
  # get plot details
  x.axis <- get.x.year.ticks(genes.mat = genes.mat)
  
  n.colors <- 100
  color.fun <- colorRampPalette(colors = c("magenta2", "grey","white"), bias = 1.5)
  col.key <- data.table("color" = color.fun(n = n.colors),
                        "value" = seq(from = min(genes.mat, na.rm = T), to = max(genes.mat, na.rm = T), along.with = 1:n.colors))
  
  # export plot and paired gene key
  
  if(nrow(genes.mat) > 0){
    pdf(file = file.path(plot.folder,"2012-timeline",paste0(taxon," ",genome,".pdf")), width = 6.5, height = 4)
    make.timeline.plot(genes.mat = genes.mat, col.key = col.key, x.axis = x.axis)
    dev.off()
  }
  
  fwrite(file = file.path(plot.folder,"2012-timeline",paste0(taxon," ",genome,".csv")), x = genes.key, sep = ",", quote = TRUE)
  
}


# ---- subset to seasons ----

genes.group <- genes[consistent.in == "Ice.On" | consistent.in == "Spring" | consistent.in == "Clearwater" | consistent.in == "Early.Summer" | consistent.in == "Late.Summer" | consistent.in == "Fall"]

if(nrow(genes.group) > 0){
  # sort genes based on occurrence patterns over time (this time need to deal with duplicates from consistent in >1 group)
  genes.wide <- dcast(data = genes.group, formula = gene + consistent.in + ko + ko.description + module + module.description + path + pathway.description ~ sample + season, value.var = "positive.selection.pvalue")
  genes.wide <- genes.wide[order(consistent.in)]
  genes.key <- genes.wide[ ,.(consistent.in = paste(consistent.in, collapse = "; ")), by = .(gene,ko,ko.description,module,module.description,path,pathway.description)]
  genes.mat <- as.matrix(x = genes.wide[ ,-c(2:8)], rownames = 1)
  genes.mat <- genes.mat[!duplicated(row.names(genes.mat)), , drop = F]
  if(!all.equal(row.names(genes.mat), genes.key$gene)){cat("\nTHERE'S A BUG IN LINE 316- key order != matrix order\n")}
  
  # Too complicated to order by season chunk, because so many in multiple seasons
  # order all together, like with overall genes
  # Hopefully this will cluster like seasons together anyway
  if (nrow(genes.mat) >= 2){
    index <- order.genes.by.appearance.pattern(genes.mat = genes.mat)
    genes.mat <- genes.mat[index, , drop = F]
    genes.key <- genes.key[index]
  }
  
  # go back from p-val = 1 to high p-vals are NA
  genes.mat[genes.mat == 1] <- NA
  
  # make seasons plot ----
  
  # sort dates
  index <- order.samples.by.date(genes.mat = genes.mat)
  genes.mat <- genes.mat[ ,index, drop = F]
  
  index <- order.samples.by.season(genes.mat = genes.mat)
  genes.mat <- genes.mat[ ,index, drop = F]
  
  # get plot details
  x.axis <- get.x.season.ticks(genes.mat = genes.mat)
  
  n.colors <- 100
  color.fun <- colorRampPalette(colors = c("magenta2", "grey","white"), bias = 1.5)
  col.key <- data.table("color" = color.fun(n = n.colors),
                        "value" = seq(from = min(genes.mat, na.rm = T), to = max(genes.mat, na.rm = T), along.with = 1:n.colors))
  
  # export plot and paired gene key
  
  if(nrow(genes.mat) > 0){
    pdf(file = file.path(plot.folder,"seasons-seasons",paste0(taxon," ",genome,".pdf")), width = 6.5, height = 4)
    make.season.plot(genes.mat = genes.mat, col.key = col.key, x.axis = x.axis)
    dev.off() 
  }
  
  fwrite(file = file.path(plot.folder,"seasons-seasons",paste0(taxon," ",genome,".csv")), x = genes.key, sep = ",", quote = TRUE)
  
}


