# RRR

# summarize the spread of matching vs lagging organisms

# ggplot sucks for barplots, can't accept wide data

# ---- set up ----

options(scipen = 9999)

library(data.table)
library(lubridate)
library(ggplot2)
library(patchwork)
library(ggnewscale)

s.cor <- fread(file = "data/2023-04-04_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-SEASONAL.tsv")
a.cor <- fread(file = "data/2023-04-04_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-ANNUAL.tsv")

tax <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
tax <- as.data.table(tax)

rel.abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
colnames(rel.abund)[3] <- "abund.perc"

cor.cutoff <- .5

# ---- parse season data ----

rel.abund <- rel.abund[ ,.(sum.abund = sum(abund.perc)), by = Genome]
tax <- merge(x = rel.abund, y = tax, by.x = "Genome", by.y = "bin.full.name", all = F)
tax[ ,phylum := sub(pattern = "p__", replacement = "", x = phylum)]
tax[is.na(phylum) | phylum == "",phylum := "Unclassified"]

s.cor <- merge(x = s.cor, y = tax, by.x = "genome", by.y = "Genome", all.x = TRUE, all.y = FALSE)
a.cor <- merge(x = a.cor, y = tax, by.x = "genome", by.y = "Genome", all.x = TRUE, all.y = FALSE)

s.cor[cor >= cor.cutoff, match := "match"]
s.cor[cor <= -cor.cutoff, match := "offset"]
s.cor[cor > -cor.cutoff & cor < cor.cutoff, match := "other"]
s.cor[ , match := factor(match, levels = c("match","offset","other"))]

a.cor[cor >= cor.cutoff, match := "match"]
a.cor[cor <= -cor.cutoff, match := "offset"]
a.cor[cor > -cor.cutoff & cor < cor.cutoff, match := "other"]
a.cor[ , match := factor(match, levels = c("match","offset","other"))]

# ---- seasonal phylum breakdown ----

s.phy <- s.cor[ , .(count = .N, abund = sum(sum.abund)), by = .(phylum, match)]

count.perc.s.phy <- dcast(data = s.phy, formula = phylum ~ match, value.var = "count")
count.perc.s.phy <- as.matrix(count.perc.s.phy, rownames = T)
count.perc.s.phy[is.na(count.perc.s.phy)] <- 0
count.perc.s.phy <- count.perc.s.phy / rowSums(count.perc.s.phy) * 100
count.perc.s.phy <- as.data.table(x = count.perc.s.phy, keep.rownames = "phylum")
count.perc.s.phy <- melt(data = count.perc.s.phy, id.vars = "phylum", variable.name = "match",value.name = "count.perc")

abund.perc.s.phy <- dcast(data = s.phy, formula = phylum ~ match, value.var = "abund")
abund.perc.s.phy <- as.matrix(abund.perc.s.phy, rownames = T)
abund.perc.s.phy[is.na(abund.perc.s.phy)] <- 0
abund.perc.s.phy <- abund.perc.s.phy / rowSums(abund.perc.s.phy) * 100
abund.perc.s.phy <- as.data.table(x = abund.perc.s.phy, keep.rownames = "phylum")
abund.perc.s.phy <- melt(data = abund.perc.s.phy, id.vars = "phylum", variable.name = "match",value.name = "abund.perc")

# ---- seasonal phylum breakdown ----

a.phy <- a.cor[ , .(count = .N, abund = sum(sum.abund)), by = .(phylum, match)]

count.perc.a.phy <- dcast(data = a.phy, formula = phylum ~ match, value.var = "count")
count.perc.a.phy <- as.matrix(count.perc.a.phy, rownames = T)
count.perc.a.phy[is.na(count.perc.a.phy)] <- 0
count.perc.a.phy <- count.perc.a.phy / rowSums(count.perc.a.phy) * 100
count.perc.a.phy <- as.data.table(x = count.perc.a.phy, keep.rownames = "phylum")
count.perc.a.phy <- melt(data = count.perc.a.phy, id.vars = "phylum", variable.name = "match",value.name = "count.perc")

abund.perc.a.phy <- dcast(data = a.phy, formula = phylum ~ match, value.var = "abund")
abund.perc.a.phy <- as.matrix(abund.perc.a.phy, rownames = T)
abund.perc.a.phy[is.na(abund.perc.a.phy)] <- 0
abund.perc.a.phy <- abund.perc.a.phy / rowSums(abund.perc.a.phy) * 100
abund.perc.a.phy <- as.data.table(x = abund.perc.a.phy, keep.rownames = "phylum")
abund.perc.a.phy <- melt(data = abund.perc.a.phy, id.vars = "phylum", variable.name = "match",value.name = "abund.perc")


# ---- make plot ----

make.stupid.barplot <- function(my.column, my.data, stupid.factor.order){
  # stupid.factor.order <- s.phy[order(my.column, decreasing = TRUE), .(..my.column = sum(..my.column)), by = phylum] # couldn't get it working in the function
  stupid.factor.order <- stupid.factor.order[order(-abund)] # wtf didn't it work the first time???
  stupid.factor.order <- stupid.factor.order[ ,phylum]
  my.data <- my.data[ , phylum := factor(phylum, levels = stupid.factor.order)]
  
  p <- ggplot(data = my.data, aes(x = phylum, y = .data[[my.column]]))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, size = 12, hjust = 1),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.title = element_text(size = 12),
          # plot.background = element_rect(),
          legend.spacing = unit(x = 2, units = "cm"),
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 10))+
    scale_fill_manual(values = c("green4","gold3","grey85"),)+
    scale_y_continuous(expand = c(.01,0))+
    geom_col(aes(fill = match),color = "black", alpha = .8, width = .8, show.legend = FALSE)+
    geom_point(aes(fill = match), color = "black", shape = 22, alpha = 0)+
    guides(fill = guide_legend(title = element_blank(), byrow = TRUE, 
                               override.aes = c(alpha = .8, size = 10))) # need byrow to get legend.spacing recognized
  return(p)
}

# OH MY GOD
# since orders are fucking FACTORS and the ggplot evaluates new each time
# you can't change the order of bars differently in different panels the pull from the same fucking data
# because the data is re-evaluated when you print the ggplot

# ---- plot by seaonal phylum ----

# stick with abundance order:
stupid.factor.order <- s.phy[order(-abund), .(abund = sum(abund)), by = phylum]

p.seas.count <- make.stupid.barplot(my.column = "count", my.data = s.phy, stupid.factor.order = stupid.factor.order)

p.seas.abund <- make.stupid.barplot(my.column = "abund", my.data = s.phy, stupid.factor.order = stupid.factor.order)

p.seas.count.perc <- make.stupid.barplot(my.column = "count.perc", my.data = count.perc.s.phy, stupid.factor.order = stupid.factor.order)

p.seas.abund.perc <- make.stupid.barplot(my.column = "abund.perc", my.data = abund.perc.s.phy, stupid.factor.order = stupid.factor.order)

# ---- plot by annual phylum ----

# stick with abundance order:
stupid.factor.order <- a.phy[order(-abund), .(abund = sum(abund)), by = phylum]

p.year.count <- make.stupid.barplot(my.column = "count", my.data = a.phy, stupid.factor.order = stupid.factor.order)

p.year.abund <- make.stupid.barplot(my.column = "abund", my.data = a.phy, stupid.factor.order = stupid.factor.order)

p.year.count.perc <- make.stupid.barplot(my.column = "count.perc", my.data = count.perc.a.phy, stupid.factor.order = stupid.factor.order)

p.year.abund.perc <- make.stupid.barplot(my.column = "abund.perc", my.data = abund.perc.a.phy, stupid.factor.order = stupid.factor.order)


# ---- combine plots ----

p.seas.abund + p.seas.count + p.seas.abund.perc + p.seas.count.perc + plot_layout(nrow = 2, guides = "collect")

p.year.abund + p.year.count + p.year.abund.perc + p.year.count.perc + plot_layout(nrow = 2, guides = "collect")

# ---- save prelim choices for fig 2 ----

saveRDS(object = p.seas.abund.perc, file = "figures/2023-04-06_fig_2_ideas/seas.abund.perc.ggplot")
saveRDS(object = p.year.abund.perc, file = "figures/2023-04-06_fig_2_ideas/ann.abund.perc.ggplot")
