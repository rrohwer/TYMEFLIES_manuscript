# RRR
# manually looked at plots with correlations .3 - .9 and decided if they looked like correlations or not
# now see how my call lined up with the actual calculated numbers

library(readxl)
library(data.table)

# ---- take 1 ----

mycalls <- read_excel(path = "figures/2023-04-04_genome_div-abund_correlations/to-choose-cutoffs/2023-04-08_manually_say_match_or_not.xlsx")
thevals <- fread(file = "data/2023-04-04_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-SEASONAL.tsv")

mycalls
thevals

mycalls$genome
mycalls$genome <- sub(pattern = " \\d+ .*$", replacement = "", x = mycalls$genome)
mycalls$genome <- sub(pattern = "^.* ", replacement = "", x = mycalls$genome)

x <- merge(y = mycalls, x = thevals, by = "genome", all.y = T, all.x = F)

x$col <- "grey"
x$col[x$call != "n"] <- "red3"  

# cor isn't working, b/c data is too smothed and noisy ones look correlated
stripchart(x = x$cor, method = "jitter", bg = x$col, pch = 21, vertical = T)  
boxplot(x = list(x$cor[x$call == "n"], x$cor[x$call == "y"]), notch = T, names = c("n","y")) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$cor[x$call == "n"], x$cor[x$call == "y"]), add = T, method = "jitter", vertical = T, pch = 21) 

# adj R2 is definitely better
stripchart(x = x$adj.R2, method = "jitter", bg = x$col, pch = 21, vertical = T)  
boxplot(x = list(x$adj.R2[x$call == "n"], x$adj.R2[x$call == "y"]), notch = T, names = c("n","y")) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$adj.R2[x$call == "n"], x$adj.R2[x$call == "y"]), add = T, method = "jitter", vertical = T, pch = 21) 

# slope seems to have no bearing
stripchart(x = x$slope, method = "jitter", bg = x$col, pch = 21, vertical = T, log = "y")  
boxplot(x = list(x$slope[x$call == "n"], x$slope[x$call == "y"]), notch = T, names = c("n","y"), ylim = c(-.1,1.5)) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$adj.R2[x$call == "n"], x$adj.R2[x$call == "y"]), add = T, method = "jitter", vertical = T, pch = 21) 

plot(x = x$adj.R2, y = x$slope, col = x$col)

# ----take 2 ----
 
# seasonal, but not averaged across years this time. .3 -.9 correlations:

mycalls <- read_excel(path = "figures/2023-04-04_genome_div-abund_correlations/to-choose-cutoffs/2023-04-08_manual_choices_take2.xlsx")
thevals <- fread(file = "data/2023-04-08_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-SEASONAL.tsv")

mycalls
thevals

mycalls$genome
mycalls$genome <- sub(pattern = " \\d+ .*$", replacement = "", x = mycalls$genome)
mycalls$genome <- sub(pattern = "^.* ", replacement = "", x = mycalls$genome)

x <- merge(y = mycalls, x = thevals, by = "genome", all.y = T, all.x = F)

x$col <- "grey"
x$col[x$choice != "n"] <- "red3"  
    
# cor isn't working, b/c data is too smothed and noisy ones look correlated
stripchart(x = x$cor, method = "jitter", bg = x$col, pch = 21, vertical = T)  
boxplot(x = list(x$cor[x$choice == "n"], x$cor[x$choice == "y"]), notch = T, names = c("n","y")) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$cor[x$choice == "n"], x$cor[x$choice == "y"]), add = T, method = "jitter", vertical = T, pch = 21) 

# adj R2 is definitely better
stripchart(x = x$adj.R2, method = "jitter", bg = x$col, pch = 21, vertical = T)  
boxplot(x = list(x$adj.R2[x$choice == "n"], x$adj.R2[x$choice == "y"]), notch = T, names = c("n","y")) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$adj.R2[x$choice == "n"], x$adj.R2[x$choice == "y"]), add = T, method = "jitter", vertical = T, pch = 21) 

# slope seems to have no bearing
stripchart(x = x$slope, method = "jitter", bg = x$col, pch = 21, vertical = T, log = "y")  
boxplot(x = list(x$slope[x$choice == "n"], x$slope[x$choice == "y"]), notch = T, names = c("n","y"), ylim = c(-.1,1.5)) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$adj.R2[x$choice == "n"], x$adj.R2[x$choice == "y"]), add = T, method = "jitter", vertical = T, pch = 21) 

plot(x = x$adj.R2, y = x$slope, col = x$col)

x[choice == "n" & cor > .7]


# ----take 3 ----

# seasonal, only genomes with anova-signif season diffs

mycalls <- read_excel(path = "figures/2023-04-04_genome_div-abund_correlations/to-choose-cutoffs/2023-04-09_manual_choices_take3.xlsx")
thevals <- fread(file = "data/2023-04-09_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-SEASONAL.tsv")

mycalls
thevals

mycalls$genome
mycalls$genome <- sub(pattern = " \\d+ .*$", replacement = "", x = mycalls$genome)
mycalls$genome <- sub(pattern = "^.* ", replacement = "", x = mycalls$genome)

x <- merge(y = mycalls, x = thevals, by = "genome", all.y = T, all.x = F)

x$col <- adjustcolor("white",alpha = 0)
x$col[x$choice == "n"] <- "grey50"  
x$col[x$choice == "m"] <- "green2"  
x$col[x$choice == "o"] <- "red2"
  
# cor looks better now
par(mar = c(.1,4,2,.1))
stripchart(x = x$cor, method = "jitter", bg = x$col, pch = 21, vertical = T, jitter = .45, xlim = c(.5,1.5), axes = F, cex = 2)  
box()
axis(side = 2, at = seq.int(from = -1, to = 1, by = .2), las = 2, cex.axis = 2)
mtext(text = "Pearson correlation", side = 3, cex = 2)

abline(h = .3)
abline(h = -.3)
boxplot(x = list(x$cor[x$choice == "n"], x$cor[x$choice == "m"], x$cor[x$choice == "o"]), notch = T, names = c("n","m","o")) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$cor[x$choice == "n"], x$cor[x$choice == "m"], x$cor[x$choice == "o"]), add = T, method = "jitter", vertical = T, pch = 21) 

# adj R2 is definitely better
stripchart(x = x$adj.R2, method = "jitter", bg = x$col, pch = 21, vertical = T, jitter = .45, xlim = c(.5,1.5))  
boxplot(x = list(x$adj.R2[x$choice == "n"], x$adj.R2[x$choice == "m"], x$adj.R2[x$choice == "o"]), notch = T, names = c("n","m","o")) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$adj.R2[x$choice == "n"], x$adj.R2[x$choice == "m"], x$adj.R2[x$choice == "o"]), add = T, method = "jitter", vertical = T, pch = 21) 

# slope seems to have no bearing
stripchart(x = x$slope, method = "jitter", bg = x$col, pch = 21, vertical = T, jitter = .45, xlim = c(.5,1.5))  
boxplot(x = list(x$slope[x$choice == "n"], x$slope[x$choice == "m"], x$slope[x$choice == "o"]), notch = T, names = c("n","m","j"), ylim = c(-.1,1.5)) #If the notches of two plots do not overlap this is ‘strong evidence’ that the two medians differ
stripchart(x = list(x$slope[x$choice == "n"], x$slope[x$choice == "m"], x$slope[x$choice == "o"]), add = T, method = "jitter", vertical = T, pch = 21) 

par(mar = c(4.5,5.5,.2,.5))
plot(x = x$adj.R2, y = x$slope, bg = adjustcolor(x$col,.8), pch = 21, cex = abs(x$cor * 4), axes = F, ann = F)
axis(side = 1, at = seq.int(from = 0, to = .8, by = .1), cex.axis = 2)
axis(side = 2, at = seq.int(from = -1, to = .3, by = .1), las = 2, cex.axis = 2)
mtext(text = "Adjusted R2", side = 1, cex = 2, line = 3)
mtext(text = "Slope", side = 2, line = 4, cex = 2)
text(x = .6, y = -.8, label = "size ~ correlation", cex = 2)
box()

abline(v = .25) # want R2 > .25

par(mar = c(4.5,5.5,.2,.5))
plot(x = x$adj.R2, y = abs(x$cor), bg = adjustcolor(x$col,.8), pch = 21, axes = F, ann = F)
axis(side = 1, at = seq.int(from = -.1, to = .8, by = .1), cex.axis = 2)
axis(side = 2, at = seq.int(from = 0, to = 1, by = .1), las = 2, cex.axis = 2)
mtext(text = "Adjusted R2", side = 1, line = 3, cex = 2)
mtext(text = "abs(Pearson Correlation)", side = 2, line = 4, cex = 2)
box()

abline(v = .25)
abline(h = .3)

x[choice == "n" & cor > .7]
x[choice == "n" & cor < -.7]

y <- copy(x)
y[ , cor := round(cor, digits = 1), ]
y <- y[ , .N, by = .(cor, choice)]
y[choice == "n", N]
y <- y[order(cor)]

par(mfrow = c(3,1))
barplot(height = y[choice == "n", N], names.arg = y[choice == "n", cor], col = "grey50")
barplot(height = y[choice == "m", N], names.arg = y[choice == "m", cor], col = "green2")
barplot(height = y[choice == "o", N], names.arg = y[choice == "o", cor], col = "red2")

