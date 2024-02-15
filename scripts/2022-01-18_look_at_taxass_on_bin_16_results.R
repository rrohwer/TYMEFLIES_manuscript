
tax <- read.csv(file = "data/2022-01-05_IMG_IDs_for_Actino_bins/2022-01-14_taxass/data/ActinoBins16S.97.80.80.taxonomy")
tax <- as.matrix(tax)

pull.p <- function(tax){
  pvals <- sub(pattern = "^.*\\(", replacement = "", x = tax)
  pvals <- sub(pattern = "\\)$", replacement = "", x = pvals)
  pvals <- as.numeric(pvals)
  return(pvals)
}
pvals <- apply(tax, 2, pull.p)

pull.tax <- function(tax){
  tax <- sub(pattern = "\\(.*\\)$", replacement = "", x = tax)
  return(tax)
}
tax <- apply(X = tax, MARGIN = 2, FUN = pull.tax)
tax[1:5, ]
tax <- tax[ ,c(2:8,1)]
colnames(tax)

alphabetical <- order(tax[,1], tax[,2], tax[,3], tax[,4], tax[,5], tax[,6], tax[,7])
tax <- tax[alphabetical, ]

otus <- rep.int(x = 1, times = nrow(tax))
bins <- list(NULL)
bins$Kingdom <- aggregate(x = otus, by = list(tax[ ,1]), FUN = sum)
bins$Phylum <- aggregate(x = otus, by = list(tax[ ,2]), FUN = sum)
bins$Class <- aggregate(x = otus, by = list(tax[ ,3]), FUN = sum)
bins$Order <- aggregate(x = otus, by = list(tax[ ,4]), FUN = sum)
bins$Lineage <- aggregate(x = otus, by = list(tax[ ,5]), FUN = sum)
bins$Clade <- aggregate(x = otus, by = list(tax[ ,6]), FUN = sum)
bins$Tribe <- aggregate(x = otus, by = list(tax[ ,7]), FUN = sum)


par(mar = c(5,4,1,0))
bar.spots <- barplot(height = bins$Phylum$x, col = "cornflowerblue")
text(x = bar.spots + .25, y = 0, labels = paste(bins$Phylum$Group.1, "  "), srt = 30, adj = 1, xpd = T)
mtext(text = "TaxAss Phylum", side = 3, cex = 1.2)
mtext(text = "(started with all IMG Actinos with a 16S)", side = 3, line = -1)
mtext(text = "Num Bins", side = 2, line = 3)

par(mar = c(5,4,1,0))
bar.spots <- barplot(height = bins$Class$x, col = "cornflowerblue")
text(x = bar.spots + .25, y = 0, labels = paste(bins$Class$Group.1, "  "), srt = 30, adj = 1, xpd = T)
mtext(text = "TaxAss Class", side = 3, cex = 1.2)
mtext(text = "(started with all IMG Actinos with a 16S)", side = 3, line = -1)
mtext(text = "Num Bins", side = 2, line = 3)

par(mar = c(5,4,1,0))
text.col <- rep("black", length(bins$Order$Group.1))
text.col[5] <- "cornflowerblue"
bar.spots <- barplot(height = bins$Order$x, col = "cornflowerblue")
text(x = bar.spots + .35, y = 0, labels = paste(sub("unclassified","U",bins$Order$Group.1), "  "), srt = 30, adj = 1, xpd = T, cex = .8,
     col = text.col)
text(x = bar.spots[5], y = bins$Order$x[5] + 5, labels = bins$Order$x[5], col = "cornflowerblue")
mtext(text = "TaxAss Order (Silva + FreshTrain)", side = 3, cex = 1.2)
mtext(text = "(started with all\nIMG Actinos with a 16S)", side = 3, line = -5, at = 20)
mtext(text = "Num Bins", side = 2, line = 3)

par(mar = c(5,4,1,0))
text.col <- rep("black", length(bins$Lineage$Group.1))
text.col[2] <- "cornflowerblue"
bar.spots <- barplot(height = bins$Lineage$x, col = "cornflowerblue")
text(x = bar.spots + .5, y = 0, labels = paste(sub("unclassified","U",bins$Lineage$Group.1), "  "), 
     srt = 40, adj = 1, xpd = T, cex = .8, col = text.col)
text(x = bar.spots[2], y = bins$Lineage$x[2] + 2, labels = bins$Lineage$x[2], col = "cornflowerblue")
mtext(text = "TaxAss Lineage (Silva + FreshTrain)", side = 3, cex = 1.2)
mtext(text = "(started with all\nIMG Actinos with a 16S)", side = 3, line = -3)
mtext(text = "Num Bins", side = 2, line = 3)

par(mar = c(5.75,4,1,0))
text.col <- rep("black", length(bins$Clade$Group.1))
text.col[2:4] <- "cornflowerblue"
bar.spots <- barplot(height = bins$Clade$x, col = "cornflowerblue")
text(x = bar.spots + .5, y = 0, labels = paste(sub("unclassified","U",bins$Clade$Group.1), "  "), 
     srt = 40, adj = 1, xpd = T, cex = .7, col = text.col)
text(x = bar.spots[2:4], y = bins$Clade$x[2:4] + 2, labels = bins$Clade$x[2:4], col = "cornflowerblue")
mtext(text = "TaxAss Clade", side = 3, cex = 1.2)
mtext(text = "(started with all\nIMG Actinos with a 16S)", side = 3, line = -3)
mtext(text = "Num Bins", side = 2, line = 3)

par(mar = c(5.75,4,1,0))
text.col <- rep("black", length(bins$Tribe$Group.1))
text.col[1:8] <- "cornflowerblue"
bar.spots <- barplot(height = bins$Tribe$x, col = "cornflowerblue")
text(x = bar.spots + .5, y = 0, labels = paste(sub("unclassified","U",bins$Tribe$Group.1), "  "), 
     srt = 40, adj = 1, xpd = T, cex = .7, col = text.col)
text(x = bar.spots[1:8], y = bins$Tribe$x[1:8] + 2, labels = bins$Tribe$x[1:8], col = "cornflowerblue")
mtext(text = "TaxAss Tribe", side = 3, cex = 1.2)
mtext(text = "(started with all\nIMG Actinos with a 16S)", side = 3, at = 15,line = -3)
mtext(text = "Num Bins", side = 2, line = 3)

par(mar = c(4,4,2,0))
acI.index <- c(1:8,19,41)
bins$acI <- bins$Tribe[acI.index, ]
text.col <- rep("black", length(bins$acI$Group.1))
text.col[1:9] <- "cornflowerblue"
bar.spots <- barplot(height = bins$acI$x, col = "cornflowerblue")
text(x = bar.spots + .25, y = 0, labels = paste(sub("unclassified","U",bins$acI$Group.1), "  "), 
     srt = 40, adj = 1, xpd = T, cex = 1, col = text.col)
text(x = bar.spots, y = bins$acI$x + .25, labels = bins$acI$x, col = text.col, xpd = T)
mtext(text = "acI Tribes", side = 3, cex = 1.2)
mtext(text = "Num Bins", side = 2, line = 3)
