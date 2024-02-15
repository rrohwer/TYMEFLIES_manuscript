# RRR

plot(x = 1:10, y = 1:10, type = "n", ann = F, axes = F)

points(x = 1.5, y = 9.8, pch = 1, xpd = NA, cex = .7)
text(x = 3.25, y = 9, labels = "Amino\nAcid\nRelated", adj = 0, xpd = NA)

points(x = 1.5, y = 6.8, pch = 2, xpd = NA, cex = .7)
text(x = 3.25, y = 6, labels = "Nucleic\nAcid\nRelated", adj = 0, xpd = NA)
