# RRR
# What's the spread of LT change over different taxa?
# plot by taxon mean abundances

by.genome <- readRDS(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/LT_change_barplot-by_genome_for_outlines.rds")
by.change <- readRDS(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/LT_change_barplot-by_change_for_colors.rds")


max(colSums(by.change))

bar.spots <- barplot(by.change, col = c("orange2","magenta3","cadetblue2"), beside = F, las = 2, border = NA, names.arg = rep("",ncol(by.change)), axes = F, legend.text = F)
barplot(by.genome, col = adjustcolor("white",0), add = T, beside = F, las = 2, names.arg = rep("",ncol(by.genome)), axes = F)
text(x = bar.spots, y = -.1, labels = colnames(by.change), srt = 45, xpd = NA, adj = 1)
axis(side = 2, at = 0:5, labels = F, lwd = 1, tck = -.02, line = -.2)
axis(side = 2, at = 0:5, labels = T, lwd = 0, line = -.5, las = 2)

mtext(text = "Mean Genome Abundance (%)", side = 2, line = 1.75)

# highlight example ones

index <- which(rownames(by.genome) == "ME2015-07-03_3300042555_group6_bin161") # p__Actinobacteriota_c__Actinomycetia_o__Nanopelagicales_f__Nanopelagicaceae_g__Planktophila_s__ ME2015-07-03_3300042555_group6_bin161
y.loc <- cumsum(by.genome[ ,1])[index] - (cumsum(by.genome[ ,1])[index] - cumsum(by.genome[ ,1])[index -1]) / 2
segments(x0 = bar.spots[1] + .25, x1 = 1.7, y0 = y.loc, y1 = y.loc)
text(x = 1.75, y = y.loc, labels = "(c)\nPlanktophila", adj = 0)

index <- which(rownames(by.genome) == "ME2011-09-21_3300043464_group3_bin69") # p__Actinobacteriota_c__Actinomycetia_o__Nanopelagicales_f__Nanopelagicaceae_g__Nanopelagicus_s__Nanopelagicus sp000294575 ME2011-09-21_3300043464_group3_bin69
y.loc <- cumsum(by.genome[ ,1])[index] - (cumsum(by.genome[ ,1])[index] - cumsum(by.genome[ ,1])[index -1]) / 2
segments(x0 = bar.spots[1] + .25, x1 = 1.7, y0 = y.loc, y1 = y.loc)
text(x = 1.75, y = y.loc, labels = "(b)\nNanopelagicus", adj = 0)

index <- which(rownames(by.genome) == "ME2005-06-22_3300042363_group2_bin84") # p__Actinobacteriota_c__Actinomycetia_o__Nanopelagicales_f__AcAMD-5_g__ATZT02_s__ATZT02 sp005789325 ME2005-06-22_3300042363_group2_bin84
y.loc <- cumsum(by.genome[ ,1])[index] - (cumsum(by.genome[ ,1])[index] - cumsum(by.genome[ ,1])[index -1]) / 2
segments(x0 = bar.spots[1] + .25, x1 = 1.7, y0 = y.loc, y1 = y.loc)
text(x = 1.75, y = y.loc, labels = "(a)\nAcAMD-5", adj = 0)

# legend

y.locs <- seq(3.75,2.5,along.with = 1:4)
rect.height <- .1
x.left <- 3.75
x.right <- 4.25
x.text <- 4.45

rect(xleft = x.left, xright = x.right, ybottom = y.locs[-1] - rect.height, ytop = y.locs[-1] + rect.height, col = c("cadetblue2","magenta3","orange2"))
text(x = x.text, y = y.locs[-1], labels = c("Gradual Change","Step Change","Disturbance/Resilience"), adj = 0)

