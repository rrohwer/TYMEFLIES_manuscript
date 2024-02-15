# RRR
# panel together Fig 2- seasonality

# half page figure is: width = 7.20472, height = 4.52756

# 7 pt font for figure text, set with pdf(..., pointsize = 7)
# get 8 pt font if default pointsize is 7 using cex = 1.142856 (this is for panel letters)

far.left <- 0
far.right <- 1
letter.width <- .01
bar.width <- .35
cor.label.width <- 0

decay.width <- far.right - bar.width - letter.width - letter.width - letter.width # minus extra letter width for more space btwn panels
corr.width <- (decay.width - letter.width - cor.label.width) / 2

left.topbar <- far.left + letter.width
right.topbar <- left.topbar + bar.width
left.decay <- right.topbar + letter.width + letter.width # plus extra letter width for more space btwn panels
right.decay <- left.decay + decay.width
left.cor1 <- far.left + letter.width + cor.label.width
right.cor1 <- left.cor1 + corr.width
left.cor2 <- right.cor1 + letter.width
right.cor2 <- left.cor2 + corr.width 
left.bar <- right.cor2 + letter.width + letter.width # plus extra letter width for more space btwn panels
right.bar <- left.bar + bar.width
left.rect.1 <- left.cor1 - letter.width
right.rect.1 <- right.cor1
left.rect.2 <- left.cor2 - letter.width
right.rect.2 <- right.cor2

all.equal(right.bar, right.decay) # check

far.top <- 1
far.bottom <- 0
letter.height <- .015
cor.label.height <- .035

bar.height <- (far.top - far.bottom - letter.height * 3) / 2 # minus extra letter width for more space btwn panels
corr.height <- (bar.height - cor.label.height) / 2

top.decay <- far.top - letter.height
bottom.decay <- top.decay - bar.height
top.abund <- bottom.decay - letter.height - letter.height # plus extra letter width for more space btwn panels
bottom.abund <- top.abund - corr.height
top.div <- bottom.abund
bottom.div <- top.div - corr.height
top.rect <- top.abund + letter.height 
bottom.rect <- bottom.div - cor.label.height

pdf(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2.pdf", width = 7.20472, height = 4.52756, title = "Rohwer_Figure_2", pointsize = 7, colormodel = "srgb")

# ---- draw the shaded boxes ----
par(fig = c(left.rect.1, right.rect.1, bottom.rect, top.rect))
par(mar = c(.25,.25,.25,.25))
plot(1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(xleft = 0, ybottom = 0, 10, 10, col = adjustcolor("gold2",.3), border = NA)
mtext(text = substitute(bold("c")), side = 3, line = -21, at = left.cor1 + .005, cex = 1.142856, xpd = NA, outer = T)

par(fig = c(left.rect.2, right.rect.2, bottom.rect, top.rect), new = T)
par(mar = c(.25,.25,.25,.25))
plot(1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(xleft = 0, ybottom = 0, 10, 10, col = adjustcolor("green4",.2), border = NA)
mtext(text = substitute(bold("d")), side = 3, line = -21, at = left.cor2 + .005, cex = 1.142856, xpd = NA, outer = T)

# ---- time decay ----
par(fig = c(left.decay, right.decay, bottom.decay, top.decay), new = T)
par(mar = c(3,3,0,.5))
source(file = "scripts-figures/2023-12-29_seasonal_figure/time_decay.R")
# box(which = "figure")
mtext(text = substitute(bold("b")), side = 3, line = -1.2, at = left.decay, cex = 1.142856, xpd = NA, outer = T)

# ---- seasonality barplot ----
par(fig = c(left.topbar, right.topbar, bottom.decay, top.decay), new = T)
par(mar = c(4.25,3.5,2.75,0))
source(file = "scripts-figures/2023-12-29_seasonal_figure/seasonal_proportion_barplot.R")
# box(which = "figure")
mtext(text = substitute(bold("a")), side = 3, line = -1.2, at = left.topbar, cex = 1.142856, xpd = NA, outer = T)

# ---- less diverse blooms ----
par(fig = c(left.cor1, right.cor1, bottom.abund, top.abund), new = T)
par(mar = c(.5,3.75,.2,.5))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-29_seasonal_figure/abund_panel-less_diverse_bloom.R")
# box(which = "figure")

par(fig = c(left.cor1, right.cor1, bottom.div, top.div), new = T)
par(mar = c(.5,3.75,.2,.5))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-29_seasonal_figure/diversity_panel-less_diverse_bloom.R")
# box(which = "figure")

# ---- more diverse blooms ----
par(fig = c(left.cor2, right.cor2, bottom.abund, top.abund), new = T)
par(mar = c(.5,3.75,.2,.5))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-29_seasonal_figure/abund_panel-more_diverse_bloom.R")
# box(which = "figure")

par(fig = c(left.cor2, right.cor2, bottom.div, top.div), new = T)
par(mar = c(.5,3.75,.2,.5))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-29_seasonal_figure/diversity_panel-more_diverse_bloom.R")
# box(which = "figure")

# ---- bloom diversity barplot ----
par(fig = c(left.bar, right.bar, bottom.rect, top.abund), new = T)
par(mar = c(4.25,3.5,2.75,0))
source(file = "scripts-figures/2023-12-29_seasonal_figure/bloom_diversity_barplot.R")
# box(which = "figure")
mtext(text = substitute(bold("e")), side = 3, line = -21, at = left.bar, cex = 1.142856, xpd = NA, outer = T)


dev.off()
