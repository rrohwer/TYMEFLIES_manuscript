# RRR
# panel together Fig 5- acI-B deeper dive

# STOP!! DON'T ADJUST PANELS WITHOUT THINKING HARDER ON THIS ONE! the NMDS aspect ratio needs to be maintained.

# half page figure is: width = 7.20472, height = 4.52756

# 7 pt font for figure text, set with pdf(..., pointsize = 7)
# get 8 pt font if default pointsize is 7 using cex = 1.142856 (this is for panel letters)

# Need to make NMDS be exactly square! well, actually exactly twice as wide as tall
# fig.width * nmds.width = nmds.inches
7.20472 * 0.485 # = 3.494289
# nmds.inches / fig.height = nmds.height
3.494289 / 4.52756 # = 0.7717819 this gives a square
0.7717819 / 2 # 0.385891 <- set height to this this maintains aspect ratio of 50:100

far.left <- 0
far.right <- 1
letter.width <- .01
nmds.width <- 0.485
heat.legendwidth <- .068

timeline.width <- far.right - nmds.width - letter.width - letter.width - letter.width - letter.width # minus extra 2 letter widths for more space btwn panels
heat.width <- nmds.width - heat.legendwidth

left.nmds <- far.left + letter.width
right.nmds <- left.nmds + nmds.width
left.timeline <- right.nmds + letter.width + letter.width + letter.width # plus extra 2 letter widths for more space btwn panels
right.timeline <- left.timeline + timeline.width
left.heat <- far.left + letter.width 
right.heat <- left.heat + heat.width 
left.heatlegend <- right.heat
right.heatlegend <- left.heatlegend + heat.legendwidth 

all.equal(right.heatlegend, right.nmds) # check

far.top <- 1
far.bottom <- 0
letter.height <- .015
nmds.height <- 0.385891
timeaxis.height <- .05
nmds.legendheight <- .05

timeline.height <- (far.top - far.bottom - timeaxis.height - letter.height * 4) / 4 
heat.height <- far.top - far.bottom - nmds.height - nmds.legendheight - letter.height * 3 # extra letter height for more spacing 

top.abund <- far.top - letter.height
bottom.abund <- top.abund - timeline.height
top.div <- bottom.abund - letter.height
bottom.div <- top.div - timeline.height
top.newsnv <- bottom.div - letter.height
bottom.newsnv <- top.newsnv - timeline.height
top.sel <- bottom.newsnv - letter.height
bottom.sel <-  top.sel - timeline.height
top.axis <- bottom.sel
bottom.axis <- far.bottom
top.nmds <- far.top - letter.height
bottom.nmds <- top.nmds - nmds.height
top.nmdslegend <- bottom.nmds 
bottom.nmdslegend <- top.nmdslegend - nmds.legendheight
top.heat <- bottom.nmdslegend - letter.height - letter.height # extra spacing
bottom.heat <- far.bottom

all.equal(far.bottom, bottom.sel - timeaxis.height) # check
all.equal(far.bottom, top.heat - heat.height) # check


pdf(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5.pdf", width = 7.20472, height = 4.52756, title = "Rohwer_Figure_5", pointsize = 7, colormodel = "srgb")

# ---- abund ----

par(fig = c(left.timeline, right.timeline, bottom.abund, top.abund))
par(mar = c(.25,4,0,.5))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/abund_timeline.R")
# box(which = "figure")
mtext(text = substitute(bold("b")), side = 3, line = -1.25, at = left.timeline - .004, cex = 1.142856, xpd = NA, outer = T)

# ---- div ----

par(fig = c(left.timeline, right.timeline, bottom.div, top.div), new = T)
par(mar = c(.25,4,0,.5))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/div_timeline.R")
# box(which = "figure")
mtext(text = substitute(bold("c")), side = 3, line = -10.75, at = left.timeline - .004, cex = 1.142856, xpd = NA, outer = T)

# ---- new SNVs ----

par(fig = c(left.timeline, right.timeline, bottom.newsnv, top.newsnv), new = T)
par(mar = c(.25,4,0,.5))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/newSNVs_timeline.R")
# box(which = "figure")
mtext(text = substitute(bold("d")), side = 3, line = -19.75, at = left.timeline - .004, cex = 1.142856, xpd = NA, outer = T)

# ---- sel ----

par(fig = c(left.timeline, right.timeline, bottom.sel, top.sel), new = T)
par(mar = c(.25,4,0,.5))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/selection_timeline.R")
# box(which = "figure")
mtext(text = substitute(bold("e")), side = 3, line = -29, at = left.timeline - .004, cex = 1.142856, xpd = NA, outer = T)

# ---- nmds ----

par(fig = c(left.nmds, right.nmds, bottom.nmds, top.nmds), new = T)
par(mar = c(1.25,1.25,.1,.1))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/nmds_plot.R")
# box(which = "figure")
mtext(text = substitute(bold("a")), side = 3, line = -1.25, at = left.nmds, cex = 1.142856, xpd = NA, outer = T)

# ---- nmds legend ----

par(fig = c(left.nmds, right.nmds, bottom.nmdslegend, top.nmdslegend), new = T)
par(mar = c(0,0,0,0))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/nmds_legend.R")
# box(which = "figure")

# ---- heatmap ----

par(fig = c(left.heat, right.heat, bottom.heat, top.heat), new = T)
par(mar = c(3,3.75,0,1))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/consistent_selection_heatmap.R")
# box(which = "figure")
mtext(text = substitute(bold("f")), side = 3, line = -19.25, at = left.nmds, cex = 1.142856, xpd = NA, outer = T)

# ---- heatmap legend (scale) ----

par(fig = c(left.heatlegend, right.heatlegend, bottom.heat, top.heat), new = T)
par(mar = c(13,.5,3.75,3))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/heatmap_legend.R")
# box(which = "figure")

# ---- heatmap legend (asterisks) ----

par(fig = c(left.heatlegend, right.heatlegend, bottom.heat, top.heat), new = T)
par(mar = c(0,.5,9,.5))
source(file = "scripts-figures/2024-01-02_acI-B_evolution_figure/heatmap_legend_asterisks.R")
# box(which = "figure")

# ----

dev.off()
