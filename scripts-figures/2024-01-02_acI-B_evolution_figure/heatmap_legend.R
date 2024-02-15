# RRR
# for fig 5f

library(data.table)

genes <- readRDS(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/consistently_selected_genes-ME2011-09-21_3300043464_group3_bin69.rds")
key <- fread(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/consistently_selected_genes_KEY-ME2011-09-21_3300043464_group3_bin69.csv")

n.colors <- 100
color.fun <- colorRampPalette(colors = c("magenta2", "grey","white"), bias = 1.5)
col.key <- data.table("color" = color.fun(n = n.colors),
                      "value" = seq(from = min(genes, na.rm = T), to = max(genes, na.rm = T), along.with = 1:n.colors))


plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 1)
axis(side = 4, at = c(.05,.15,.25), xpd = T, lwd = 0, las = 2, line = -.5)
axis(side = 4, at = c(.05,.15,.25), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
mtext(text = "Positive\nSelection\np-value", adj = 0, line = 1)
