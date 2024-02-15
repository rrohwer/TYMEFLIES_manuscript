# RRR
# new SNVs plot for ESA, no season colors, base R version 

library(data.table)
library(lubridate)

acI.A <- readRDS("figures/2023-07-18_tracking_new_SNVs/acI_A.rds")
acI.B <- readRDS("figures/2023-07-18_tracking_new_SNVs/acI_B.rds")
acI.C <- readRDS("figures/2023-07-18_tracking_new_SNVs/acI_C.rds")

info <- acI.C

# ---- parse data ----
info <- info[order(Date)]

# ---- plot set-up ----

add.magic.boxes <- function(){
  box(which = "inner", col="red", lwd = 3)
  box(which = "outer", col="blue", lwd = 3)
  box(which = "plot", col="purple", lwd = 3)
  box(which = "figure", col="orange", lwd = 3)
  # inner margins (orange to purple) are mar/plt/mai
  # outer margins (red to blue) are oma/omd/omi
  par.settings <- par(no.readonly = TRUE)
  cat("change distance btwn ORANGE and PURPLE with mar  (mar = ", par.settings$mar, ")\n")
  cat("change distance btwn BLUE and RED with oma       (oma = ", par.settings$oma, ")\n")
}

fill.under.lines <- function(X, Y, YAxisMin, Color, xpd = F){
  poly.x <- c(min(X), X, max(X))
  poly.y <- c(YAxisMin, Y, YAxisMin )
  polygon(x = poly.x, y = poly.y, col = Color, border = NA, xpd = xpd)
}

x.ticks <- parse_date_time(x = 2000:2019, orders = "y")

# acI-B abund ----
# pdf(file = "figures/2023-08-05_ESA_plots/new_SNVs_acI-B_long-term.pdf",width = 6.5, height = 3)
# pdf(file = "figures/2023-08-05_ESA_plots/new_SNVs_acI-A_long-term.pdf",width = 6.5, height = 3)
pdf(file = "figures/2023-08-05_ESA_plots/new_SNVs_acI-C_long-term.pdf",width = 6.5, height = 3)
par(mar = c(1.75,3.5,1,0))
plot(x = info$Date, y = info$New, type = "n", ann = F, axes = F)
fill.under.lines(X = info$Date, Y = info$New, YAxisMin = 0, Color = adjustcolor(col = "#8e99a7", alpha.f = .5))
lines(x = info$Date, y = info$New, col = "#8e99a7")
points(x = info$Date, y = info$New, pch = 16, col = adjustcolor("#8e99a7", alpha.f = .7))
axis(side = 1, at = c(min(x.ticks), max(info$Date)), lwd = 1, lwd.ticks = 0, labels = F)
axis(side = 1, at = x.ticks[seq.int(1,20,2)], labels = F, tck = -.025, lwd = 0, lwd.ticks = 1)
axis(side = 1, at = x.ticks[seq.int(1,20,2)], labels = year(x.ticks)[seq.int(1,20,2)], lwd = 0, las = 1, line = -.5)
axis(side = 2, labels = F, line = -.75, tck = -.025)
axis(side = 2, lwd = F, las = 2, line = -1.25)
mtext(text = "New SNVs (count)", side = 2, line = 2.5)
dev.off()
# ----
