# RRR

library(tidyr)
library(lubridate)

plot.folder <- "figures/2023-02-16_div_vs_cov/"

genome.info <- readRDS("data/2023-12-12_instrain_results_TIMEOUT/genome_info_combined.rds")
bin.info <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")
bin.info <- bin.info[bin.info$winner, ]

colnames(genome.info)

# make some matrices

nuc.div <- pivot_wider(data = genome.info, id_cols = "genome", names_from = "sample", values_from = "nucl_diversity") |>
  as.data.frame()
row.names(nuc.div) <- nuc.div[ ,1]
nuc.div <- nuc.div[ ,-1] |>
  as.matrix()
nuc.div.max <- apply(X = nuc.div, MARGIN = 1, FUN = max, na.rm = T)
nuc.div.perc <- nuc.div / nuc.div.max * 100

cov.med <- pivot_wider(data = genome.info, id_cols = "genome", names_from = "sample", values_from = "coverage_median") |>
  as.data.frame()
row.names(cov.med) <- cov.med[ ,1]
cov.med <- cov.med[ ,-1] |>
  as.matrix()
cov.med.max <- apply(X = cov.med, MARGIN = 1, FUN = max, na.rm = T)
cov.med.perc <- cov.med / cov.med.max * 100

index <- which(cov.med.max >= 100) # 70 is top 200 bins
length(index)
top.cov.med <- cov.med[index, ]
top.cov.med.perc <- cov.med.perc[index, ]
top.nuc.div <- nuc.div[index, ]
top.nuc.div.perc <- nuc.div.perc[index, ]

top.orgs <- data.frame("genome" = row.names(top.nuc.div.perc), "row.order" = 1:nrow(top.nuc.div.perc))
top.orgs <- merge(x = top.orgs, y = bin.info, by.x = "genome", by.y = "bin.full.name", all.x = T, all.y = F)
top.orgs <- top.orgs[order(top.orgs$row.order), ]
all.equal(top.orgs$genome, row.names(top.nuc.div.perc))



# plot everything
pdf(file = file.path(plot.folder,"div_vs_cov.pdf"), width = 4, height = 4)
par(mar = c(4,4,1,1))
plot(x = c(min(nuc.div, na.rm = T), max(nuc.div, na.rm = T)), y = c(min(cov.med, na.rm = T), max(cov.med, na.rm = T)), type = "n", ann = F)
mtext(text = "Nucleotide Diversity", side = 1, line = 2.5)
mtext(text = "Median Coverage", side = 2, line = 2.5)
for(g in 1:nrow(nuc.div)){
  points(x = nuc.div[g, ], y = cov.med[g, ], col = adjustcolor("darkblue",.2))
}
dev.off()

pdf(file = file.path(plot.folder,"div_vs_cov-perc_max.pdf"), width = 4, height = 4)
par(mar = c(4,4,1,1))
plot(x = c(min(nuc.div.perc, na.rm = T), max(nuc.div.perc, na.rm = T)), y = c(min(cov.med.perc, na.rm = T), max(cov.med.perc, na.rm = T)), type = "n", ann = F)
mtext(text = "Nucleotide Diversity (% of max)", side = 1, line = 2.5)
mtext(text = "Median Coverage (% of max)", side = 2, line = 2.5)
for(g in 1:nrow(nuc.div)){
  points(x = nuc.div.perc[g, ], y = cov.med.perc[g, ], col = adjustcolor("darkblue",.2))
}
dev.off()

# plot top 200 only

pdf(file = file.path(plot.folder,"div_vs_cov-top_orgs.pdf"), width = 4, height = 4)
par(mar = c(4,4,1,1))
plot(x = c(min(top.nuc.div, na.rm = T), max(top.nuc.div, na.rm = T)), y = c(min(top.cov.med, na.rm = T), max(top.cov.med, na.rm = T)), type = "n", ann = F)
mtext(text = "Nucleotide Diversity", side = 1, line = 2.5)
mtext(text = "Median Coverage", side = 2, line = 2.5)
for(g in 1:nrow(top.nuc.div)){
  points(x = top.nuc.div[g, ], y = top.cov.med[g, ], col = adjustcolor("darkblue",.2))
}
dev.off()

pdf(file = file.path(plot.folder,"div_vs_cov-perc_max-top_orgs.pdf"), width = 4, height = 4)
par(mar = c(4,4,1,1))
plot(x = c(min(top.nuc.div.perc, na.rm = T), max(top.nuc.div.perc, na.rm = T)), y = c(min(top.cov.med.perc, na.rm = T), max(top.cov.med.perc, na.rm = T)), type = "n", ann = F)
mtext(text = "Nucleotide Diversity (% of max)", side = 1, line = 2.5)
mtext(text = "Median Coverage (% of max)", side = 2, line = 2.5)
for(g in 1:nrow(top.nuc.div)){
  points(x = top.nuc.div.perc[g, ], y = top.cov.med.perc[g, ], col = adjustcolor("darkblue",.2))
}
dev.off()

# color by taxa

for (p in unique(top.orgs$phylum)){
  pdf(file = file.path(plot.folder,"phylum",paste0(p,".pdf")), width = 4, height = 4)
  par(mar = c(4,4,2,1))
  plot(x = c(min(top.nuc.div.perc, na.rm = T), max(top.nuc.div.perc, na.rm = T)), y = c(min(top.cov.med.perc, na.rm = T), max(top.cov.med.perc, na.rm = T)), type = "n", ann = F)
  mtext(text = "Nucleotide Diversity (% of max)", side = 1, line = 2.5)
  mtext(text = "Median Coverage (% of max)", side = 2, line = 2.5)
  mtext(text = p, side = 3, line = .5)
  col.key <- rep.int(adjustcolor("grey40",.2),times = nrow(top.nuc.div.perc))
  col.key[top.orgs$phylum == p] <- "red4"
  for(g in 1:nrow(top.nuc.div)){
     points(x = top.nuc.div.perc[g, ], y = top.cov.med.perc[g, ], col = adjustcolor(col.key[g],.5))
  }
  dev.off()
}

# actinos and proteos and cyanos and planctomycetes were pretty crowded in that view

for (p in unique(top.orgs$phylum)){
  for (f in unique(top.orgs$family[top.orgs$phylum == p])){
    pdf(file = file.path(plot.folder,"family",paste0(p,"_",f,".pdf")), width = 4, height = 4)
    par(mar = c(4,4,2,1))
    plot(x = c(min(top.nuc.div.perc, na.rm = T), max(top.nuc.div.perc, na.rm = T)), y = c(min(top.cov.med.perc, na.rm = T), max(top.cov.med.perc, na.rm = T)), type = "n", ann = F)
    mtext(text = "Nucleotide Diversity (% of max)", side = 1, line = 2.5)
    mtext(text = "Median Coverage (% of max)", side = 2, line = 2.5)
    mtext(text = paste(p,f,sep = ", "), side = 3, line = .5)
    col.key <- rep.int(adjustcolor("grey40",.2),times = nrow(top.nuc.div.perc))
    col.key[top.orgs$family == f] <- "red4"
    pch.key <- rep(1,length(col.key))
    pch.key[col.key == "red4"] <- 21
    for(g in 1:nrow(top.nuc.div)){
      points(x = top.nuc.div.perc[g, ], y = top.cov.med.perc[g, ], col = col.key[g], pch = pch.key[g], bg = col.key[g])
    }
    dev.off()
  }
}


# look at top orgs individually

for(g in 1:nrow(top.nuc.div.perc)){
  taxonomy <- paste(bin.info[bin.info$bin.full.name == row.names(top.nuc.div)[g],c("phylum","class","order","family","genus","species","num.in.cluster")], collapse = " ")
  pdf(file = file.path(plot.folder,"top_orgs",paste0(taxonomy," - ",g,".pdf")), width = 4, height = 4)
  par(mar = c(4,4,3,1))
  plot(x = top.nuc.div.perc[g, ], y = top.cov.med.perc[g, ], col = adjustcolor("darkblue",.9), ann = F)
  mtext(text = "Nucleotide Diversity (% of max)", side = 1, line = 2.5)
  mtext(text = "Median Coverage (% of max)", side = 2, line = 2.5)
  mtext(text = row.names(top.nuc.div.perc)[g], side = 3, line = .5)
  mtext(text = taxonomy, cex = .7, side = 3, line = 1.5)
  dev.off()
}

