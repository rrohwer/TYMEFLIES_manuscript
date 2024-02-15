# RRR
# plot how many reads map to the bins, broken down by season

mapped <- readRDS("figures/2023-12-24_paper_figures/Rohwer_Figure_1_data/percent_mapped.rds")




min(unlist(mapped))
max(unlist(mapped))

stripchart(x = mapped, vertical = T, method = "jitter", jitter = .25, axes = F, ann = F, pch = 16, col = adjustcolor("grey40",.5), cex = .7, add = F)
boxplot(mapped, pch = NA, ann = F, axes = F, lwd = 1, lty = 1, add = T, 
        col = c(adjustcolor("snow3",.3),
                adjustcolor("tan4",.6),
                adjustcolor("cornflowerblue",.6),
                adjustcolor("chartreuse4",.6),
                adjustcolor("purple",.4),
                adjustcolor("hotpink2",.5)), 
        # border = c(adjustcolor("snow3",1),
        #            adjustcolor("tan4",1),
        #            adjustcolor("cornflowerblue",1),
        #            adjustcolor("chartreuse4",1),
        #            adjustcolor("purple",1),
        #            adjustcolor("hotpink2",1))
        border = "black")
box()
axis(side = 2, at = seq(20,65,15), labels = F, lwd = 0, lwd.ticks = 1, tck = -.025)
axis(side = 2, at = seq(20,65,15), labels = T, lwd = 0, las = 2, line = -.5)
mtext(text = "Reads Mapped (%)", side = 2, line = 1.9)
