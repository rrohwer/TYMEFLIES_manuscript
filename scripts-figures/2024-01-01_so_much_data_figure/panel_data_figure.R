# RRR
# panel together Fig 1- ALL the samples, ALL the bins!

# half page figure is: width = 7.20472, height = 4.52756
# quarter page figure is: width = 3.50394, height = 4.52756
# half page vertical figure would be: width = 3.50394, height = 9.724409

# 7 pt font for figure text, set with pdf(..., pointsize = 7)
# get 8 pt font if default pointsize is 7 using cex = 1.142856 (this is for panel letters)

far.left <- 0
far.right <- 1
letter.width <- .015
qual.width <- .2
rank.names.width <- .22

mapped.width <- far.right - qual.width * 2 - letter.width * 5 # minus extra letter width for more space btwn panels
rank.width <- (far.right - rank.names.width - letter.width * 2) / 2

left.dates <- far.left + letter.width
right.dates <- far.right
left.legend <- left.dates
right.legend <- right.dates
left.qual.1 <- left.dates
right.qual.1 <- left.qual.1 + qual.width 
left.qual.2 <- right.qual.1 + letter.width
right.qual.2 <- left.qual.2 + qual.width
left.mapped <- right.qual.2 + letter.width + letter.width + letter.width # plus extra letter width for more space btwn panels
right.mapped <- left.mapped + mapped.width
left.names <- far.left + letter.width
right.names <- left.names + rank.names.width
left.rank.1 <- right.names
right.rank.1 <- left.rank.1 + rank.width
left.rank.2 <- right.rank.1 + letter.width
right.rank.2 <- left.rank.2 + rank.width

all.equal(right.mapped, far.right)
all.equal(right.rank.2, far.right)

far.top <- 1
far.bottom <- 0
letter.height <- .015
dates.height <- .33
legend.height <- .05
# qual.height <- .33

remaining.height <- far.top - far.bottom - dates.height - legend.height - letter.height * 6 # minus extra letter height for more space btwn panels
mapped.height <- remaining.height * .5 
bar.height <- remaining.height * .5
qual.height <- mapped.height

top.dates <- far.top - letter.height
bottom.dates <- top.dates - dates.height
top.legend <- bottom.dates - letter.height
bottom.legend <- top.legend - legend.height
top.qual <- bottom.legend - letter.height - letter.height # extra space
bottom.qual <- top.qual - qual.height
top.mapped <- top.qual
bottom.mapped <- top.mapped - mapped.height
top.bar <- bottom.mapped - letter.height - letter.height # extra space
bottom.bar <- far.bottom

all.equal(far.bottom, top.bar - bar.height)

pdf(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_1.pdf", width = 3.50394, height = 4.52756, title = "Rohwer_Figure_1", pointsize = 7, colormodel = "srgb")

# sample dates ----

par(fig = c(left.dates, right.dates, bottom.dates, top.dates))
par(mar = c(1.25,2.5,0,.2))
source(file = "scripts-figures/2024-01-01_so_much_data_figure/plot_sample_dates.R")
# box(which = "figure")
mtext(text = substitute(bold("a")), side = 3, line = -1, at = left.dates, cex = 1.142856, xpd = NA, outer = T)

# season legend ----

par(fig = c(left.legend, right.legend, bottom.legend, top.legend), new = T)
par(mar = c(.1,.1,.1,.1))
source(file = "scripts-figures/2024-01-01_so_much_data_figure/make_season_legend.R")
# box(which = "figure")

# completeness ----

par(fig = c(left.qual.1, right.qual.1, bottom.qual, top.qual), new = T)
par(mar = c(1.4,2.6,.1,.1))
source(file = "scripts-figures/2024-01-01_so_much_data_figure/plot_bin_completeness.R")
# box(which = "figure")
mtext(text = substitute(bold("b")), side = 3, line = -17.5, at = left.qual.1, cex = 1.142856, xpd = NA, outer = T)

# contamination ----

par(fig = c(left.qual.2, right.qual.2, bottom.qual, top.qual), new = T)
par(mar = c(1.4,2.6,.1,.1))
source(file = "scripts-figures/2024-01-01_so_much_data_figure/plot_bin_contamination.R")
# box(which = "figure")

# mapped reads ----

par(fig = c(left.mapped, right.mapped, bottom.mapped, top.mapped), new = T)
par(mar = c(.1,2.8,.1,.2))
source(file = "scripts-figures/2024-01-01_so_much_data_figure/plot_reads_mapped.R")
# box(which = "figure")
mtext(text = substitute(bold("c")), side = 3, line = -17.5, at = left.mapped, cex = 1.142856, xpd = NA, outer = T)

# rank abund 16S ----

par(fig = c(left.rank.1, right.rank.1, bottom.bar, top.bar), new = T)
par(mar = c(.1,.1,2.25,.7))
source(file = "scripts-figures/2024-01-01_so_much_data_figure/plot_16S_rank_abund.R")
# box(which = "figure")
mtext(text = substitute(bold("d")), side = 3, line = -29, at = left.rank.1, cex = 1.142856, xpd = NA, outer = T)

# rank abund MAGs ----

par(fig = c(left.rank.2, right.rank.2, bottom.bar, top.bar), new = T)
par(mar = c(.1,.1,2.25,.7))
source(file = "scripts-figures/2024-01-01_so_much_data_figure/plot_MAGs_rank_abund.R")
# box(which = "figure")
mtext(text = substitute(bold("e")), side = 3, line = -29, at = left.rank.2, cex = 1.142856, xpd = NA, outer = T)

dev.off()
