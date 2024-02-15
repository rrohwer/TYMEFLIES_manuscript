# RRR
# I realized I can skip correlating the clusters, and instead just directly correlate the *genomes*
# was doing clustering because I thought I'd match them by hand, but that was too hard to do.
# not sure if worth it to pursue analyzing the clusters anymore.
# steal from that script to group and trim the dataset.
# go back and trin less harshly- don't want to exclude all cyanobacteria!

# ---- set up ----

library(data.table)
library(lubridate)
library(compositions) # for CLR transform
library(ggplot2)
library(patchwork)
library(limony)
library(ggnewscale)

genome.info <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")

bin.info <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds") 
bin.info <- as.data.table(bin.info)

rel.abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
colnames(rel.abund)[3] <- "abund.perc"
rel.abund <- rel.abund[Genome != "unmapped", .(sample = Sample, genome = Genome, abund.perc)]

data("seasons")

tax <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

output.plot.folder <- "figures/2023-04-04_genome_div-abund_correlations/"
output.data.folder <- "data/2023-04-04_genome_abund-div_correlations/"

breadth.cutoff <- .5
season.cutoff <- 3 # must be in this many seasons to be included in seasonal analysis
season.cutoff.2 <- 5 # must be in a season in this many years to calculate an average value for that season
annual.cutoff <- 10 # must be in this many years (in any season) to be included in the annual analysis

# ---- make abundance zero based on breadth ----

colnames(genome.info)
hist(genome.info$breadth / genome.info$breadth_expected)

x <- genome.info[ , .N, by = genome] 
nrow(x) # 2855 genomes
sum(x$N) # present 255892 total times

genome.info <- genome.info[ (breadth / breadth_expected) >= breadth.cutoff, ] # **** is this too lenient?? ****

x <- genome.info[ , .N, by = genome] 
nrow(x) # 2853 genomes
sum(x$N) # present 240278 total times
255892 - 240278 # 15,000 occurrences got removed 

rel.abund <- merge(x = genome.info, y = rel.abund, by = c("genome","sample"), all.x = TRUE, all.y = FALSE)

# ---- average div and abund data to account for uneven sampling ----

# average for each season in each year
div.by.seas <- genome.info[order(genome, year, season), .(nuc.div = mean(nucl_diversity, na.rm = T)), by = .(genome, year, season)]
abund.by.seas <- rel.abund[order(genome, year, season), .(abund.perc = mean(abund.perc, na.rm = T)), by = .(genome, year, season)]

# average seasons across all years- only include if at least season.cutoff.2 years have observations of that season
div.by.seas.overlay <- div.by.seas[ , .(.N, nuc.div = mean(nuc.div, na.rm = T)), by = .(genome, season)]
abund.by.seas.overlay <- abund.by.seas[ , .(.N, abund.perc = mean(abund.perc, na.rm = T)), by = .(genome, season)]

div.by.seas.overlay <- div.by.seas.overlay[N >= season.cutoff.2]
abund.by.seas.overlay <- abund.by.seas.overlay[N >= season.cutoff.2]
div.by.seas.overlay[ , N := NULL]
abund.by.seas.overlay[ , N := NULL]

# average years across all seasons
div.by.year.avgs <- div.by.seas[ , .(nuc.div = mean(nuc.div, na.rm = T)), by = .(genome, year)]
abund.by.year.avgs <- abund.by.seas[ , .(abund.perc = mean(abund.perc, na.rm = T)), by = .(genome, year)]

# make into matrix
d.seasonal <- dcast(data = div.by.seas.overlay, formula = genome ~ season, value.var = "nuc.div")
d.seasonal <- as.matrix(d.seasonal, rownames = T)

a.seasonal <- dcast(data = abund.by.seas.overlay, formula = genome ~ season, value.var = "abund.perc")
a.seasonal <- as.matrix(a.seasonal, rownames = T)

d.annual <- dcast(data = div.by.year.avgs, formula = genome ~ year, value.var = "nuc.div")
d.annual <- as.matrix(d.annual, rownames = T)

a.annual <- dcast(data = abund.by.year.avgs, formula = genome ~ year, value.var = "abund.perc")
a.annual <- as.matrix(a.annual, rownames = T)

# ---- filter based on pres/abs ----

# include only genomes present in at least this many seasons / year 
index <- apply(!is.na(d.seasonal), 1, sum)
index <- index >= season.cutoff
d.seasonal <- d.seasonal[index, ]
a.seasonal <- a.seasonal[index, ]

# include only genomes present in at least this many years
index <- which(colnames(d.annual) == "2019") # exclude 2019 partial year
d.annual <- d.annual[, -index]
a.annual <- a.annual[, -index]

index <- apply(!is.na(d.annual), 1, sum)
index <- index >= annual.cutoff
d.annual <- d.annual[index, ]
a.annual <- a.annual[index, ]

# ---- make data normally distributed for pearson correlation ----

# pearson corr wants normal data
# why does dividing by the mean make div more normal??
# CLR includes a log transform, makes sense the abund is lognormal so it also makes it normal
# CLR is better than just a log, because technically shouldn't do correlations on compositions

# normalize by diversity mean values
my.means <- apply(d.seasonal, 1, mean, na.rm = T)
d.seasonal <- d.seasonal / my.means
hist(d.seasonal, breaks = 100)

my.means <- apply(d.annual, 1, mean, na.rm = T)
d.annual <- d.annual / my.means
hist(d.annual, breaks = 100)

# do CLR transform to get abundance more normal

# check rows vs columns this package expects:
x <- matrix(data = c(2,20,5,50,7,70,3,30,1,10,2,20), nrow = 4)
rownames(x) <- letters[1:4]
colnames(x) <- LETTERS[1:3]
x # each row is a composition vector (rows sum to 1)
clr(x) # the rows with the same composition have the same CLR values

# in my data, each column is a composition vector (unclosed)
a.seasonal <- a.seasonal |>
  t() |>
  clr() |>
  t() |>
  as.matrix()
hist(a.seasonal)

a.annual <- a.annual |>
  t() |>
  clr() |>
  t() |>
  as.matrix()
hist(a.annual)  

# ---- calc correlation for each genome ----

d.seasonal <- t(d.seasonal)
d.annual <- t(d.annual)
a.seasonal <- t(a.seasonal)
a.annual <- t(a.annual)

cor.seasonal <- stats::cor(x = d.seasonal, y = a.seasonal, method = "pearson", use = "pairwise.complete.obs") # I think it's faster to do the whole matrix than loop through only genome:same genome comparisons.
cor.annual <- stats::cor(x = d.annual, y = a.annual, method = "pearson", use = "pairwise.complete.obs")

cor.seasonal <- as.data.table(x = cor.seasonal, keep.rownames = "genome1")
cor.annual <- as.data.table(x = cor.annual, keep.rownames = "genome1")

cor.seasonal <- melt(data = cor.seasonal, id.vars = 1, variable.name = "genome2", value.name = "cor")
cor.annual <- melt(data = cor.annual, id.vars = 1, variable.name = "genome2", value.name = "cor")

cor.seasonal <- cor.seasonal[genome1 == genome2]
cor.annual <- cor.annual[genome1 == genome2]

hist(cor.seasonal$cor, breaks = 100, xlim = c(-1,1)) # pretty even
hist(cor.annual$cor, breaks = 100, xlim = c(-1,1)) # kind of normal

# ---- what about a linear model to determine correlation? ----

# pete suggested a model, which has the benefit of a p-value

lm.seasonal <- data.table("genome" = colnames(a.seasonal),"adj.R2" = 0,"slope" = 0)
all.equal(colnames(a.seasonal), colnames(d.seasonal))
for (g in 1:ncol(d.seasonal)){
  my.lm <- lm(d.seasonal[ ,g] ~ a.seasonal[ ,g])
  my.lm <- summary(my.lm)
  lm.seasonal[g, adj.R2 := my.lm$adj.r.squared]
  lm.seasonal[g, slope := my.lm$coefficients[2,1]]
}

lm.annual <- data.table("genome" = colnames(a.annual),"adj.R2" = 0,"slope" = 0)
all.equal(colnames(a.annual), colnames(d.annual))
for (g in 1:ncol(d.annual)){
  my.lm <- lm(d.annual[ ,g] ~ a.annual[ ,g])
  my.lm <- summary(my.lm)
  lm.annual[g, adj.R2 := my.lm$adj.r.squared]
  lm.annual[g, slope := my.lm$coefficients[2,1]]
}

# combine all the metrics into 1 table

cor.seasonal <- merge(x = cor.seasonal, y = lm.seasonal, by.x = "genome1", by.y = "genome")
cor.seasonal[ ,genome2 := NULL]

cor.annual <- merge(x = cor.annual, y = lm.annual, by.x = "genome1", by.y = "genome")
cor.annual[ ,genome2 := NULL]

# ---- plot seasonal correlations - choose cutoff for "signif" correlations ----

# get x-axis coords for average season boundaries
data("seasons") # from limony package
seasons[ ,-1] <- apply(seasons[ ,-1], MARGIN = 2, FUN = yday)
seasons$Ice.On[seasons$Ice.On < 100] <- seasons$Ice.On[seasons$Ice.On < 100] + 365
season.avs <- colMeans(seasons[-1])
season.midpoints <- season.avs 
season.midpoints[1] <- (season.midpoints[2] - season.midpoints[1]) / 2 + season.midpoints[1]
season.midpoints[2] <- (season.midpoints[3] - season.midpoints[2]) / 2 + season.midpoints[2]
season.midpoints[3] <- (season.midpoints[4] - season.midpoints[3]) / 2 + season.midpoints[3]
season.midpoints[4] <- (season.midpoints[5] - season.midpoints[4]) / 2 + season.midpoints[4]
season.midpoints[5] <- (season.midpoints[6] - season.midpoints[5]) / 2 + season.midpoints[5]
season.midpoints[6] <- (season.midpoints[1] - 1) / 2 + 1
season.midpoints <- data.table("season" = names(season.midpoints),"midpoint" = season.midpoints, "start" = season.avs)

# add into single tables for ggplot
div.by.seas.overlay <- merge(x = season.midpoints, y = div.by.seas.overlay, by = "season")
abund.by.seas.overlay <- merge(x = season.midpoints, y = abund.by.seas.overlay, by = "season")

# make factor for ordered ggplot color assignments
genome.info[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
rel.abund[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
div.by.seas.overlay[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
abund.by.seas.overlay[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]

# get tax for plot labels
tax$full.taxonomy <- paste(tax$phylum,tax$class,tax$order,tax$family,tax$genus,tax$species)
tax$bin.label <- paste0(tax$bin.full.name, " (",tax$num.in.cluster,")")
tax$subtitle <- paste(tax$full.taxonomy, tax$bin.label, sep = "\n")

# get month labels for x axis
x.ax.ticks <- yday(parse_date_time(x = 1:12, orders = "m"))
x.ax.labs <- month(parse_date_time(x = 1:12, orders = "m"), label = T, abbr = T)

for (g in cor.seasonal$genome1){
  p.seas.div <- ggplot(data = genome.info[genome == g], aes(x = yday, y = nucl_diversity))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_point(aes(color = season), shape = 16, size = 3, alpha = .9, show.legend = F)+
    geom_line(data = div.by.seas.overlay[genome == g], aes(x = midpoint, y = nuc.div), alpha = 1, lwd = 3)+
    geom_point(data = div.by.seas.overlay[genome == g], aes(x = midpoint, y = nuc.div, col = season), shape = 15, size = 10, alpha = .6)+
    guides(color = guide_legend(title = "Season", override.aes = c(alpha = 1)))+
    scale_x_continuous(name = element_blank(), breaks = x.ax.ticks, labels = x.ax.labs)+
    ylab("Nucleotide Diversity (pi)")
  
  p.seas.abund <- ggplot(data = rel.abund[genome == g], aes(x = yday, y = abund.perc))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_point(aes(color = season), shape = 16, size = 3, alpha = .9, show.legend = F)+
    geom_line(data = abund.by.seas.overlay[genome == g], aes(x = midpoint, y = abund.perc), alpha = 1, lwd = 3)+
    geom_point(data = abund.by.seas.overlay[genome == g], aes(x = midpoint, y = abund.perc, col = season), shape = 15, size = 10, alpha = .6)+
    guides(color = guide_legend(title = "Season", override.aes = c(alpha = 1)))+
    scale_x_continuous(name = element_blank(), breaks = x.ax.ticks, labels = x.ax.labs)+
    ylab("Relatice Abundance (CLR-transformed %)")
  
  # plots without stats for choosing cutoffs
  p.seas <- p.seas.abund + p.seas.div + plot_layout(nrow = 2, guides = "collect") +
    plot_annotation(title = g,
                    theme = theme(plot.title = element_text(hjust = .5, size = 20)))

  if (cor.seasonal[genome1 == g, cor] >= .3 & cor.seasonal[genome1 == g, cor] <= .9){
    pdf(file = file.path(output.plot.folder,"to-choose-cutoffs/seasonal",paste(tax$full.taxonomy[tax$bin.full.name == g], tax$bin.full.name[tax$bin.full.name == g], tax$num.in.cluster[tax$bin.full.name == g], ".pdf")),
        width = 10, height = 9)
    print(p.seas)
    dev.off()
  }
  
  # # plots with all the data and stats
  # p.seas <- p.seas.abund + p.seas.div + plot_layout(nrow = 2, guides = "collect") + 
  #   plot_annotation(title = paste("COR =", round(x = cor.seasonal[genome1 == g, cor], digits = 2), 
  #                                 "     adj R2 =", round(x = cor.seasonal[genome1 == g, adj.R2], digits = 2), 
  #                                 " slope =", round(x = cor.seasonal[genome1 == g, slope], digits = 2)), 
  #                   subtitle =  tax$subtitle[tax$bin.full.name == g], 
  #                   theme = theme(plot.title = element_text(hjust = .5, size = 20)))
  # 
  # pdf(file = file.path(output.plot.folder,"seasonal",paste(tax$full.taxonomy[tax$bin.full.name == g], tax$bin.full.name[tax$bin.full.name == g], tax$num.in.cluster[tax$bin.full.name == g], ".pdf")),
  #     width = 10, height = 9)
  # print(p.seas)
  # dev.off()
}

# ---- plot annual correlations - choose cutoff for "signif" correlations ----

# make factor for ggplot coloring
genome.info[ ,invasion := factor(invasion, levels = c("none","spiny","zebra"))]

div.by.year.avgs$invasion <- "none"
div.by.year.avgs[year >= 2010, invasion := "spiny"]
div.by.year.avgs[year >= 2016, invasion := "zebra"]
div.by.year.avgs[ ,invasion := factor(invasion, levels = c("none","spiny","zebra"))]

abund.by.year.avgs$invasion <- "none"
abund.by.year.avgs[year >= 2010, invasion := "spiny"]
abund.by.year.avgs[year >= 2016, invasion := "zebra"]
abund.by.year.avgs[ ,invasion := factor(invasion, levels = c("none","spiny","zebra"))]

year.midpoints <- data.table("year.midpoint" = parse_date_time(x = paste(2000:2018, "183"), orders = "yj"), "year" = 2000:2018)
div.by.year.avgs <- merge(x = year.midpoints, y = div.by.year.avgs, by = "year")
abund.by.year.avgs <- merge(x = year.midpoints, y = abund.by.year.avgs, by = "year")

# make date formats match to share plot axis
genome.info[ ,date := parse_date_time(x = date, orders = "ymd")]
rel.abund[ ,date := parse_date_time(x = date, orders = "ymd")]

for (g in cor.annual$genome1){
  p.time.div <- ggplot(data = genome.info[genome == g], aes(x = date, y = nucl_diversity))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_point(aes(color = season), shape = 16, size = 1.5, alpha = .9)+
    new_scale_color()+
    scale_color_manual(values = c("steelblue","orange2","red3"))+
    geom_line(data = div.by.year.avgs[genome == g], aes(x = year.midpoint, y = nuc.div), lwd = 3)+
    geom_point(data = div.by.year.avgs[genome == g], aes(x = year.midpoint, y = nuc.div, color = invasion), shape = 15, size = 10, alpha = .6)+
    scale_x_continuous(name = element_blank(), breaks = parse_date_time(x = 2000:2019, orders = "y"), labels = 2000:2019)+
    ylab("Nucleotide Diversity (pi)")
  
  p.time.abund <- ggplot(data = rel.abund[genome == g], aes(x = date, y = abund.perc))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_point(aes(color = season), shape = 16, size = 1.5, alpha = .9)+
    new_scale_color()+
    scale_color_manual(values = c("steelblue","orange2","red3"))+
    geom_line(data = abund.by.year.avgs[genome == g], aes(x = year.midpoint, y = abund.perc), lwd = 3)+
    geom_point(data = abund.by.year.avgs[genome == g], aes(x = year.midpoint, y = abund.perc, color = invasion), shape = 15, size = 10, alpha = .6)+
    scale_x_continuous(name = element_blank(), breaks = parse_date_time(x = 2000:2019, orders = "y"), labels = 2000:2019)+
    ylab("Relatice Abundance (CLR-transformed %)")
  
  # blided plots to choose the cutoff
  p.time <- p.time.abund + p.time.div + plot_layout(nrow = 2, guides = "collect") +
    plot_annotation(title = g,
                    theme = theme(plot.title = element_text(hjust = .5, size = 20)))

  pdf(file = file.path(output.plot.folder,"to-choose-cutoffs/annual",paste(tax$full.taxonomy[tax$bin.full.name == g], tax$bin.full.name[tax$bin.full.name == g], tax$num.in.cluster[tax$bin.full.name == g], ".pdf")),
      width = 10, height = 9)
  print(p.time)
  dev.off()
  
  # # pretty plots with all the info
  # p.time <- p.time.abund + p.time.div + plot_layout(nrow = 2, guides = "collect") + 
  #   plot_annotation(title = paste("COR =", round(x = cor.annual[genome1 == g, cor], digits = 2)),
  #                   subtitle =  tax$subtitle[tax$bin.full.name == g], 
  #                   theme = theme(plot.title = element_text(hjust = .5, size = 20)))
  # 
  # pdf(file = file.path(output.plot.folder,"long-term",paste(tax$full.taxonomy[tax$bin.full.name == g], tax$bin.full.name[tax$bin.full.name == g], tax$num.in.cluster[tax$bin.full.name == g], ".pdf")),
  #     width = 10, height = 9)
  # print(p.time)
  # dev.off()
  
}




# ---- save data ----

saveRDS(object = div.by.seas.overlay, file = file.path(output.data.folder,"diversity_grouped_by_season.rds"))
saveRDS(object = abund.by.seas.overlay, file = file.path(output.data.folder, "abundance_grouped_by_season.rds"))

saveRDS(object = div.by.year.avgs, file = file.path(output.data.folder,"diversity_grouped_by_year.rds"))
saveRDS(object = abund.by.year.avgs, file = file.path(output.data.folder, "abundance_grouped_by_year.rds"))

colnames(cor.seasonal)[1] <- "genome"
colnames(cor.annual)[1] <- "genome"

fwrite(x = cor.seasonal, file = file.path(output.data.folder,"correlations_btwn_diversity_and_abundance-SEASONAL.tsv"))
fwrite(x = cor.annual, file = file.path(output.data.folder,"correlations_btwn_diversity_and_abundance-ANNUAL.tsv"))

# ---- save plot idea for a rough fig 2 ----
# 
# g <- "ME2002−09−12pf_3300042556_group1_bin222" # acI-C? season match
# g <- "ME2005−06−22_3300042363_group2_bin36" # acI-B? season offset
# g <- "ME2011−09−04_3300044729_group3_bin142" # acI-A? season offset
# g <- "ME2011−09−21_3300043464_group3_bin69" # acI-B? neither but LT change in nuc div
# g <- "ME2015−03−07_3300042541_group6_bin53" # acI-A LT offset
# 
# g <- "ME2002-09-12pf_3300042556_group1_bin222" # acI-C? season match
# cor.seasonal[genome1 == g]
# # re-run season plots
# p.season.match <- p.seas.abund + p.seas.div + plot_layout(nrow = 2, guides = "collect") + 
#   plot_annotation(title = "acI-C", 
#                   theme = theme(plot.title = element_text(hjust = .5, size = 20)))
# 
# g <- "ME2011-09-04_3300044729_group3_bin142" # acI-A? season offset
# cor.seasonal[genome1 == g]
# # re-run season plots
# p.season.offset <- p.seas.abund + p.seas.div + plot_layout(nrow = 2, guides = "collect") + 
#   plot_annotation(title = "acI-A", 
#                   theme = theme(plot.title = element_text(hjust = .5, size = 20)))
# 
# g <- "ME2011-09-21_3300043464_group3_bin69" # acI-B? neither but LT change in nuc div
# cor.annual[genome1 == g]
# # re-run season plots
# p.LT.other <- p.time.abund + p.time.div + plot_layout(nrow = 2, guides = "collect") + 
#   plot_annotation(title = "acI-B", 
#                   theme = theme(plot.title = element_text(hjust = .5, size = 20)))
# 
# saveRDS(object = p.season.match, file = "figures/2023-04-06_fig_2_ideas/acI-C_season_match.ggplot")
# saveRDS(object = p.season.offset, file = "figures/2023-04-06_fig_2_ideas/acI-A_season_offset.ggplot")
# saveRDS(object = p.LT.other, file = "figures/2023-04-06_fig_2_ideas/acI-B_long-term_other.ggplot")
