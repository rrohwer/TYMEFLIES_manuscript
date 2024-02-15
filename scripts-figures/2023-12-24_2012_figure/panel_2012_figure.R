# RRR

# one script that souced each figure script and panels them together
# that way data is uploaded separately, so easier to troubleshoot
# and paneling is more clear and easier to adjust

set.seed(seed = 20240122)

far.left <- 0
far.right <- 1
temp.width <- .26
box.width <- .09
arrow.width <- 0 #.01
letter.width <- .01
rect.width <- .01 #.02

nano.left <- far.left + letter.width
nano.right <- far.right - (temp.width + (box.width * 3) + arrow.width + (letter.width * 4) + rect.width * 4)
temp.left <- nano.right + letter.width + rect.width
temp.right <- temp.left + temp.width
ice.left <- temp.right + letter.width
ice.right <- ice.left + box.width
zoop.left <- ice.right + letter.width
zoop.right <- zoop.left + box.width
arrow.left <- zoop.right + rect.width
arrow.right <- arrow.left + arrow.width 
bio.left <- arrow.right + rect.width + letter.width
bio.right <- bio.left + box.width

hot.rect.left <- nano.right
hot.rect.right <- zoop.right + rect.width
bio.rect.left <- arrow.right
bio.rect.right <- bio.right + rect.width


far.top <- 1
far.bottom <- 0
letter.height <- .025
rect.height.top <- .04
rect.height.bottom <- .04
temp.height <- (far.top / 2) - letter.height  - rect.height.top - rect.height.bottom

nano.top <- far.top - letter.height
nano.bottom <- far.bottom
temp.top <- far.top - letter.height - rect.height.top
temp.bottom <- temp.top - temp.height
dis.top <- temp.bottom - rect.height.top - rect.height.bottom - letter.height
dis.bottom <- dis.top - temp.height
bio.top <- temp.top - rect.height.top
bio.bottom <- bio.top - temp.height

hot.rect.bottom <- temp.bottom - rect.height.bottom
hot.rect.top <- far.top
dry.rect.bottom <- far.bottom
dry.rect.top <- dis.top + letter.height + rect.height.top

# ---- make plot ----

# half page figure is: width = 7.20472, height = 4.52756

pdf(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4.pdf", width = 7.20472, height = 3, title = "Rohwer_Figure_4", pointsize = 7, colormodel = "srgb")


# beeswarm abrupt change in Nanopelagicaceae plot ----

par(fig = c(nano.left, nano.right, nano.bottom, nano.top))
par(mar = c(8,2,0,.5))
source(file = "scripts-figures/2023-12-24_2012_figure/abrupt_snv_change.R")
# box(which = "figure")
mtext(text = substitute(bold("a")), side = 3, line = -1, at = nano.left, cex = 1.142856, xpd = NA, outer = T) # 8 pt font if default pointsize is 7

# hot rectangle section ----

par(fig = c(hot.rect.left, hot.rect.right, hot.rect.bottom, hot.rect.top), new = T)
# par(mar = c(.25,1,.25,.25))
par(mar = c(.25,.25,.25,.25))
plot(1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(xleft = 0, ybottom = 0, 10, 10, col = adjustcolor("red4",.1), border = NA)
mtext(text = "Hot Year", side = 3, line = -1, adj = .2, outer = F, col = "red4")

# dry rectangle section ----

par(fig = c(hot.rect.left, hot.rect.right, dry.rect.bottom, dry.rect.top), new = T)
# par(mar = c(.25,1,.25,.25))
par(mar = c(.25,.25,.25,.25))
plot(1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(xleft = 0, ybottom = 0, xright = 10, ytop = 10, col = adjustcolor("yellow3",.1), border = NA)
mtext(text = "Dry Year", side = 3, line = -1, adj = .2, outer = F, col = "yellow4")

# bio rectangle section ----

par(fig = c(bio.rect.left, bio.rect.right, dry.rect.bottom, hot.rect.top), new = T)
# par(mar = c(.25,1,.25,.25))
par(mar = c(.25,.25,.25,.25))
plot(1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(xleft = 0, ybottom = 0, xright = 10, ytop = 10, col = adjustcolor("green4", .1), border = NA)
mtext(text = "Low", side = 3, line = -1, adj = 0, at = 3.5, outer = F, col = "green4")
mtext(text = "Productivity", side = 3, line = -2, adj = 0, at = 3.5, outer = F, col = "green4")

# # arrow ----
# 
# par(fig = c(arrow.left, arrow.right, dry.rect.top - .1, hot.rect.bottom + .1), new = T)
# par(mar = c(0,0,0,0))
# plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
# arrows(x0 = 2, x1 = 20, y0 = 5.5, y1 = 5.5, length = .07, angle = 30, xpd = NA, lwd = 2)

# temp plot ----

par(fig = c(temp.left, temp.right, temp.bottom, temp.top), new = T)
par(mar = c(1.5,3,.1,.1))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-24_2012_figure/epi_temp.R")
mtext(text = substitute(bold("b")), side = 3, line = -2, at = temp.left, cex = 1.142856, xpd = NA, outer = T) # 8 pt font if default pointsize is 7

# ice plot ----

par(fig = c(ice.left, ice.right, temp.bottom, temp.top), new = T)
par(mar = c(.1,3.5,.1,.1))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-24_2012_figure/ice_duration.R")
# box(which = "figure")
mtext(text = substitute(bold("c")), side = 3, line = -2, at = ice.left + .005, cex = 1.142856, xpd = NA, outer = T)

# zoop plot ----

par(fig = c(zoop.left, zoop.right, temp.bottom, temp.top), new = T)
par(mar = c(.1,3.5,.1,.1))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-24_2012_figure/mean_zoops.R")
# box(which = "figure")
mtext(text = substitute(bold("d")), side = 3, line = -2, at = zoop.left + .005, cex = 1.142856, xpd = NA, outer = T)

# discharge plot ----

par(fig = c(temp.left, temp.right, dis.bottom, dis.top), new = T)
par(mar = c(1.5,3,.1,.1))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-24_2012_figure/yahara_discharge.R")
mtext(text = substitute(bold("e")), side = 3, line = -14.8, at = temp.left, cex = 1.142856, xpd = NA, outer = T)

# TP plot ----

par(fig = c(ice.left, ice.right, dis.bottom, dis.top), new = T)
par(mar = c(.1,3.5,.1,.1))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-24_2012_figure/epi_mean_TP.R")
# box(which = "figure")
mtext(text = substitute(bold("f")), side = 3, line = -14.8, at = ice.left + .005, cex = 1.142856, xpd = NA, outer = T)

# SRP plot ----

par(fig = c(zoop.left, zoop.right, dis.bottom, dis.top), new = T)
par(mar = c(.1,3.5,.1,.1))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-24_2012_figure/epi_mean_SRP.R")
# box(which = "figure")
mtext(text = substitute(bold("g")), side = 3, line = -14.8, at = zoop.left + .005, cex = 1.142856, xpd = NA, outer = T)

# biomass plot ----

par(fig = c(bio.left, bio.right, bio.bottom, bio.top), new = T)
par(mar = c(.1,2.5,.1,.1))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-24_2012_figure/mean_phyto.R")
# box(which = "figure")
mtext(text = substitute(bold("h")), side = 3, line = -2 - 1.2, at = bio.left, cex = 1.142856, xpd = NA, outer = T)

# DOC plot ----

par(fig = c(bio.left, bio.right, dis.bottom, dis.top), new = T)
par(mar = c(.1,2.5,.1,.1))
plot(1:10,1:10, type = "n", axes = F, ann = F, xaxs = "i", yaxs = "i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
par(new = T)
source(file = "scripts-figures/2023-12-24_2012_figure/epi_mean_DOC.R")
# box(which = "figure")
mtext(text = substitute(bold("i")), side = 3, line = -14.8 - .2, at = bio.left, cex = 1.142856, xpd = NA, outer = T)

dev.off()
