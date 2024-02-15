# RRR

library(ggplot2)
library(patchwork)

p.bar.seas <- readRDS("figures/2023-04-06_fig_2_ideas/seas.abund.perc.ggplot")
p.bar.ann <- readRDS("figures/2023-04-06_fig_2_ideas/ann.abund.perc.ggplot")

p.match.seas <- readRDS("figures/2023-04-06_fig_2_ideas/acI-C_season_match.ggplot")
p.offset.seas <- readRDS("figures/2023-04-06_fig_2_ideas/acI-A_season_offset.ggplot")
p.other.ann <- readRDS("figures/2023-04-06_fig_2_ideas/acI-B_long-term_other.ggplot")

layout <- "
AB
AB
AC
DC
DE
DE
"

# pdf(file = "figures/2023-04-06_fig_2_ideas/Fig_2_rough_idea-pdf", width = 4.76, height = 8)
pdf(file = "figures/2023-04-06_fig_2_ideas/Fig_2_rough_idea-pdf", width = 9.52, height = 16)
p.bar.seas + p.match.seas + p.offset.seas + p.bar.ann + p.other.ann + plot_layout(design = layout, guides = "collect", widths = c(1,2), heights = c(1,1,1))
dev.off()
