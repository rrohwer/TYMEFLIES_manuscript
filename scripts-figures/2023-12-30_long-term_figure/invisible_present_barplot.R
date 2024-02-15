# RRR
# barplot showing how LT patterns often overlay seasonal patterns
# fig 3d

lt <- readRDS("figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/invisible_present_barplot_matrix.rds")

lt[1, ] <- lt[1, ] * -1

x.lim <- c(min(lt),max(lt))

bar.spots <- barplot(lt[1, ], beside = T, horiz = T, border = T, xlim = x.lim, col = c("orange2","magenta3","cadetblue2"), ann = F, axes = F, names.arg = F, space = .75)
barplot(lt[2, ], beside = T, horiz = T, border = F, xlim = x.lim, add = T, col = c("orange2","magenta3","cadetblue2"), ann = F, axes = F, names.arg = F, space = .75)
barplot(lt[2, ], beside = T, horiz = T, border = T, col = adjustcolor("black", alpha.f = .2), xlim = x.lim, add = T, ann = F, axes = F, names.arg = F, space = .75)

bar.widths <- bar.spots[2] - bar.spots[1]

text(x = x.lim[1] - 17.5, y = bar.spots, labels = c("Disturbance/Resilience","Step Change","Gradual Change"), adj = 0, xpd = NA)
text(x = c(-1),y = bar.spots[3] + bar.widths /1.2, labels = "Seasonal", xpd = NA, adj = 1)
text(x = c(1),y = bar.spots[3] + bar.widths /1.2, labels = "Not Seasonal", xpd = NA, adj = 0)
text(x = lt[1, ] - .5, y = bar.spots, labels = abs(lt[1,]), xpd = NA, adj = 1)        
text(x = lt[2, ] + .5, y = bar.spots, labels = abs(lt[2,]), xpd = NA, adj = 0)        
text(x = 0, y = bar.spots[1] - bar.widths /1, labels = "Total Genomes: 263", xpd = NA, adj = .3)

segments(x0 = 0, x1 = 0, y0 = bar.spots[1] - .5, y1 = bar.spots[3] + 2, xpd = NA)

