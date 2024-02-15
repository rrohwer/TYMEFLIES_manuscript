# RRR
# for fig 5f

library(data.table)
library(readxl)

genes <- readRDS(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/consistently_selected_genes-ME2011-09-21_3300043464_group3_bin69.rds")
key <- read_excel(path = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/consistently_selected_genes_KEY-ME2011-09-21_3300043464_group3_bin69-plus_manual_edits.xlsx")
key <- as.data.table(key)

x.years <- substr(x = colnames(genes), start = 3, stop = 6)
x.ticks <- which(!duplicated(x.years))
x.labs <- unique(x.years)

y.categories <- key[nrow(key):1, consistent.in]
y.ticks <- c(which(!duplicated(y.categories)), length(y.categories) + 1)
y.ticks <- y.ticks - .5
y.labs <- unique(y.categories)
y.labs <- sub("^p","P",y.labs)
y.labs <- sub("^o","O",y.labs)
y.lab.locs <- key[nrow(key):1, .(row.order,consistent.in)]
y.lab.locs[ ,reverse.row.order := 1:nrow(y.lab.locs)]
y.lab.locs <- y.lab.locs[ , .(y.lab.locs = mean(reverse.row.order)), by = .(consistent.in)]
y.lab.locs <- y.lab.locs$y.lab.locs

aa.locs <- key[nrow(key):1, amino.acid.related]
aa.locs <- which(aa.locs)
na.locs <- key[nrow(key):1, nucleic.acid.related]
na.locs <- which(na.locs)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("magenta2", "grey","white"), bias = 1.5)
col.key <- data.table("color" = color.fun(n = n.colors),
                      "value" = seq(from = min(genes, na.rm = T), to = max(genes, na.rm = T), along.with = 1:n.colors))

image.mat <- t(genes)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
gene.lab <- sub(pattern = "ME2011-09-21_3300043464_group3_bin69_scaffold_", replacement = "", x = colnames(image.mat))


image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F,
      col = col.key$color)
box(which = "plot", lwd = 1)

axis(side = 1, at = x.ticks, labels = F, lwd.ticks = 1, lwd = 0, tck = -.03)
axis(side = 1, at = x.ticks, labels = x.labs, lwd = 0, las = 2, hadj = .5, line = .75)

# axis(side = 2, at = 1:nrow(genes), lwd = 0, labels = gene.lab, las = 2, line = -.75, cex.axis = .4)
axis(side = 2, at = y.ticks, labels = F, lwd.ticks = 1, lwd = 0, tck = -.05)
axis(side = 2, at = y.lab.locs[1], labels = y.labs[1], lwd = 0, xpd = NA, line = -.5)
axis(side = 2, at = y.lab.locs[2:3], labels = y.labs[2:3], lwd = 0, xpd = NA, las = 2, line = -.75)

points(x = rep(nrow(image.mat),length(aa.locs)) + 10, y = aa.locs, pch = 1, xpd = NA, cex = .5)
points(x = rep(nrow(image.mat),length(na.locs)) + 10, y = na.locs, pch = 2, xpd = NA, cex = .5)

mtext(text = "Consistently Selected Genes", side = 2, line = 3, at = y.lab.locs[1])


# stats for text:
key[ ,.N,by = consistent.in]
