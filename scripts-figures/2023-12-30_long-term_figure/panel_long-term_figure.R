# RRR
# panel together Fig 3- long-term change

# half page figure is: width = 7.20472, height = 4.52756

# 7 pt font for figure text, set with pdf(..., pointsize = 7)
# get 8 pt font if default pointsize is 7 using cex = 1.142856 (this is for panel letters)

far.left <- 0
far.right <- 1
letter.width <- .0025
bar.width <- .35
rect.width <- .005

timeline.width <- far.right - bar.width - letter.width * 4 - rect.width * 2 # extra letter width for more spacing 

rect.left <- far.left
timeline.left <- far.left + rect.width + letter.width
timeline.right <- timeline.left + timeline.width
rect.right <- timeline.right + rect.width
bar.left <- rect.right + letter.width * 3 # extra letter width for more spacing 
bar.right <- bar.left + bar.width

all.equal(bar.right, far.right) # check


far.top <- 1
far.bottom <- 0
letter.height <- .0075
rect.height.top <- .015
rect.height.bottom <- .01
sidebar.height <- .2

timeline.height <- (far.top - letter.height * 3 - rect.height.top * 3 - rect.height.bottom * 3) / 3
bar.height <- (far.top - sidebar.height - letter.height * 4)  # extra letter heights for more spacing

top.grad.rect <- far.top - letter.height
top.grad <- top.grad.rect - rect.height.top
bottom.grad <- top.grad - timeline.height
bottom.grad.rect <- bottom.grad - rect.height.bottom
top.step.rect <- bottom.grad.rect - letter.height
top.step <- top.step.rect - rect.height.top
bottom.step <- top.step - timeline.height
bottom.step.rect <- bottom.step - rect.height.bottom
top.dist.rect <- bottom.step.rect - letter.height
top.dist <- top.dist.rect - rect.height.top
bottom.dist <- top.dist - timeline.height
bottom.dist.rect <- bottom.dist - rect.height.bottom
all.equal(bottom.dist.rect, far.bottom)
bottom.dist.rect <- far.bottom # rounding errors made it negative and not a valid fig input

top.sidebar <- far.top - letter.height
bottom.sidebar <- top.sidebar - sidebar.height
top.bar <- bottom.sidebar - letter.height * 3 # extra letter height for more spacing
bottom.bar <- top.bar - bar.height

all.equal(bottom.bar, far.bottom)


pdf(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_3.pdf", width = 7.20472, height = 4.52756, title = "Rohwer_Figure_3", pointsize = 7, colormodel = "srgb")

# ---- draw shaded boxes ----
par(fig = c(rect.left, rect.right, bottom.grad.rect, top.grad.rect))
par(mar = c(.1,.1,.1,.1))
plot(1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(xleft = 0, ybottom = 0, 10, 10, col = adjustcolor("cadetblue2",.3), border = NA)
mtext(text = substitute(bold("a")), side = 3, line = -1.5, at = timeline.left, cex = 1.142856, xpd = NA, outer = T)

par(fig = c(rect.left, rect.right, bottom.step.rect, top.step.rect), new = T)
par(mar = c(.1,.1,.1,.1))
plot(1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(xleft = 0, ybottom = 0, 10, 10, col = adjustcolor("magenta3",.2), border = NA)
mtext(text = substitute(bold("b")), side = 3, line = -14.5, at = timeline.left, cex = 1.142856, xpd = NA, outer = T)

par(fig = c(rect.left, rect.right, bottom.dist.rect, top.dist.rect), new = T)
par(mar = c(.1,.1,.1,.1))
plot(1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(xleft = 0, ybottom = 0, 10, 10, col = adjustcolor("orange2",.25), border = NA)
mtext(text = substitute(bold("c")), side = 3, line = -27.5, at = timeline.left, cex = 1.142856, xpd = NA, outer = T)

# ---- draw LT change examples ----

par(fig = c(timeline.left, timeline.right, bottom.grad, top.grad), new = T)
par(mar = c(1.4,3,.1,.6))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source("scripts-figures/2023-12-30_long-term_figure/gradual_change_example_plot.R")
# box(which = "figure")

par(fig = c(timeline.left, timeline.right, bottom.step, top.step), new = T)
par(mar = c(1.4,3,.1,.6))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source("scripts-figures/2023-12-30_long-term_figure/step_change_example_plot.R")
# box(which = "figure")

par(fig = c(timeline.left, timeline.right, bottom.dist, top.dist), new = T)
par(mar = c(1.4,3,.1,.6))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source("scripts-figures/2023-12-30_long-term_figure/dist_resilience_example_plot.R")
# box(which = "figure")

# ---- sideways barplot ----

par(fig = c(bar.left, bar.right, bottom.sidebar, top.sidebar), new = T)
par(mar = c(1.5,8.25,1.5,2.4))
source("scripts-figures/2023-12-30_long-term_figure/invisible_present_barplot.R")
# box(which = "figure")
mtext(text = substitute(bold("d")), side = 3, line = -1.5, at = bar.left + .005, cex = 1.142856, xpd = NA, outer = T)

# ---- bar plot ----

par(fig = c(bar.left, bar.right, bottom.bar, top.bar), new = T)
par(mar = c(5.4,2.7,.1,0))
source("scripts-figures/2023-12-30_long-term_figure/LT_change_by_taxa_barplot.R")
# box(which = "figure")
mtext(text = substitute(bold("e")), side = 3, line = -9.5, at = bar.left + .005, cex = 1.142856, xpd = NA, outer = T)

dev.off()






