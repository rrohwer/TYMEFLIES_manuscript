# RRR
# my manual choices didn't match the correlation values very well, 
# mostly I think because I said no when it looked like there was a lot of noise,
# even when the averaged lines matched pretty well.
# so I think that means I'm averaging too much. 
# re-do that script, but this time get an average value for each season in the year, 
# but don't average the years together to have 1 value / season

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
output.data.folder <- "data/2023-04-08_genome_abund-div_correlations"

breadth.cutoff <- .7
season.cutoff <- 3 # must be in this many seasons to be included in seasonal analysis
season.cutoff.2 <- 5 # must be in a season in this many years to count as in that season
annual.cutoff <- 10 # must be in this many years to be included in the annual analysis
annual.cutoff.2 <- 3 # must be in each year at least this many times to be counted as present in that year

# ---- exclude observations where breadth is much lower than expected ----

hist(genome.info$breadth / genome.info$breadth_expected)

genome.info <- genome.info[ (breadth / breadth_expected) >= breadth.cutoff, ] 

rel.abund <- merge(x = genome.info, y = rel.abund, by = c("genome","sample"), all.x = TRUE, all.y = FALSE)

# ---- remove all instances where abundance is 0 (b/c nuc div will not exist) ----

rel.abund <- rel.abund[abund.perc > 0] # breadth cutoff already removed them?

# ---- filter for seasonal analysis  ----

# average for each season in each year
d.seasonal <- rel.abund[order(genome, year, season), .(nuc.div = mean(nucl_diversity)), by = .(genome, year, season)]
a.seasonal <- rel.abund[order(genome, year, season), .(abund.perc = mean(abund.perc)), by = .(genome, year, season)]

# be in at least 3 seasons at least 5 times to include
keep <- d.seasonal[ , .N, by = .(season, genome)]
keep[ ,N := N >= season.cutoff.2]
keep <- keep[ , .(N = sum(N)), by = genome]
keep <- keep[N >= season.cutoff]
d.seasonal <- merge(x = keep[ , .(genome)], y = d.seasonal, by = "genome", all.x = T, all.y = F)
a.seasonal <- merge(x = keep[ , .(genome)], y = a.seasonal, by = "genome", all.x = T, all.y = F)

# ---- filter for annual analysis ----

# be in a year at least twice for 10 years to include
keep <- rel.abund[ , .N, by = .(year, genome)]
keep[ , N := N >= annual.cutoff.2]
keep <- keep[ , .(N = sum(N)), by = genome]
keep <- keep[N >= annual.cutoff]
keep.annual <- merge(x = keep, y = rel.abund, by = "genome", all.x = T, all.y = F)

# average all non-zero abundances across the year
d.annual <- keep.annual[order(genome, year, season), .(nuc.div = mean(nucl_diversity)), by = .(genome, year)]
a.annual <- keep.annual[order(genome, year, season), .(abund.perc = mean(abund.perc)), by = .(genome, year)]

# ---- make into matrix ----

d.seasonal <- dcast(data = d.seasonal, formula = genome ~ season + year, value.var = "nuc.div")
d.seasonal <- as.matrix(d.seasonal, rownames = T)
a.seasonal <- dcast(data = a.seasonal, formula = genome ~ season + year, value.var = "abund.perc")
a.seasonal <- as.matrix(a.seasonal, rownames = T)

d.annual <- dcast(data = d.annual, formula = genome ~ year, value.var = "nuc.div")
d.annual <- as.matrix(d.annual, rownames = T)
a.annual <- dcast(data = a.annual, formula = genome ~ year, value.var = "abund.perc")
a.annual <- as.matrix(a.annual, rownames = T)

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

hist(cor.seasonal$cor, breaks = 100, xlim = c(-1,1)) # kind of normal
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
  lm.seasonal[g, pvalue := my.lm$coefficients[4]] # is this the p-value? I'm not sure!!!
}

lm.annual <- data.table("genome" = colnames(a.annual),"adj.R2" = 0,"slope" = 0)
all.equal(colnames(a.annual), colnames(d.annual))
for (g in 1:ncol(d.annual)){
  my.lm <- lm(d.annual[ ,g] ~ a.annual[ ,g])
  my.lm <- summary(my.lm)
  lm.annual[g, adj.R2 := my.lm$adj.r.squared]
  lm.annual[g, slope := my.lm$coefficients[2,1]]
  lm.annual[g, pvalue := my.lm$coefficients[4]] # is this the p-value? I'm not sure!!!
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

# prep formats for ggplot
rel.abund <- merge(x = season.midpoints, y = rel.abund, by = "season")
rel.abund[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]

# get tax for plot labels
tax$full.taxonomy <- paste(tax$phylum,tax$class,tax$order,tax$family,tax$genus,tax$species)
tax$bin.label <- paste0(tax$bin.full.name, " (",tax$num.in.cluster,")")
tax$subtitle <- paste(tax$full.taxonomy, tax$bin.label, sep = "\n")

# get month labels for x axis
x.ax.ticks <- yday(parse_date_time(x = 1:12, orders = "m"))
x.ax.labs <- month(parse_date_time(x = 1:12, orders = "m"), label = T, abbr = T)

for (g in cor.seasonal$genome1){
  p.seas.div <- ggplot(data = rel.abund[genome == g], aes(x = yday, y = nucl_diversity))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_point(aes(color = season), shape = 16, size = 3, alpha = .9, show.legend = T)+
    guides(color = guide_legend(title = "Season", override.aes = c(alpha = 1, shape = 15, size = 5)))+
    scale_x_continuous(name = element_blank(), breaks = x.ax.ticks, labels = x.ax.labs)+
    ylab("Nucleotide Diversity (pi)")
  
  p.seas.abund <- ggplot(data = rel.abund[genome == g], aes(x = yday, y = abund.perc))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_point(aes(color = season), shape = 16, size = 3, alpha = .9, show.legend = T)+
    guides(color = guide_legend(title = "Season", override.aes = c(alpha = 1, shape = 15, size = 5)))+
    scale_x_continuous(name = element_blank(), breaks = x.ax.ticks, labels = x.ax.labs)+
    ylab("Relatice Abundance (%)")
  
  # plots without stats for choosing cutoffs
  p.seas <- p.seas.abund + p.seas.div + plot_layout(nrow = 2, guides = "collect") +
    plot_annotation(title = g,
                    theme = theme(plot.title = element_text(hjust = .5, size = 20)))
  
  if (cor.seasonal[genome1 == g, cor] >= .3 & cor.seasonal[genome1 == g, cor] <= .9){
    pdf(file = file.path(output.plot.folder,"to-choose-cutoffs-take2/seasonal",paste(tax$full.taxonomy[tax$bin.full.name == g], tax$bin.full.name[tax$bin.full.name == g], tax$num.in.cluster[tax$bin.full.name == g], ".pdf")),
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


# make date formats match axis date format
genome.info[ ,date := parse_date_time(x = date, orders = "ymd")]
rel.abund[ ,date := parse_date_time(x = date, orders = "ymd")]

for (g in cor.annual$genome1){
  p.time.div <- ggplot(data = rel.abund[genome == g], aes(x = date, y = nucl_diversity))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_point(aes(color = season), shape = 16, size = 1.5, alpha = .9)+
    guides(color = guide_legend(title = "Season", override.aes = c(alpha = 1, shape = 15, size = 5)))+
    scale_x_continuous(name = element_blank(), breaks = parse_date_time(x = 2000:2019, orders = "y"), labels = 2000:2019)+
    ylab("Nucleotide Diversity (pi)")
  
  p.time.abund <- ggplot(data = rel.abund[genome == g], aes(x = date, y = abund.perc))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_point(aes(color = season), shape = 16, size = 1.5, alpha = .9)+
    guides(color = guide_legend(title = "Season", override.aes = c(alpha = 1, shape = 15, size = 5)))+
    scale_x_continuous(name = element_blank(), breaks = parse_date_time(x = 2000:2019, orders = "y"), labels = 2000:2019)+
    ylab("Relatice Abundance (%)")
  
  # blided plots to choose the cutoff
  p.time <- p.time.abund + p.time.div + plot_layout(nrow = 2, guides = "collect") +
    plot_annotation(title = g,
                    theme = theme(plot.title = element_text(hjust = .5, size = 20)))
  
  pdf(file = file.path(output.plot.folder,"to-choose-cutoffs-take2/annual",paste(tax$full.taxonomy[tax$bin.full.name == g], tax$bin.full.name[tax$bin.full.name == g], tax$num.in.cluster[tax$bin.full.name == g], ".pdf")),
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

saveRDS(object = d.seasonal, file = file.path(output.data.folder,"diversity_filtered_by_seasonal_abundance.rds"))
saveRDS(object = a.seasonal, file = file.path(output.data.folder, "abundance_filtered_by_seasonal_abundance.rds"))

saveRDS(object = d.annual, file = file.path(output.data.folder,"diversity_filtered_by_annual_abundance.rds"))
saveRDS(object = a.annual, file = file.path(output.data.folder, "diversity_filtered_by_annual_abundance.rds"))

colnames(cor.seasonal)[1] <- "genome"
colnames(cor.annual)[1] <- "genome"

fwrite(x = cor.seasonal, file = file.path(output.data.folder,"correlations_btwn_diversity_and_abundance-SEASONAL.tsv"))
fwrite(x = cor.annual, file = file.path(output.data.folder,"correlations_btwn_diversity_and_abundance-ANNUAL.tsv"))
