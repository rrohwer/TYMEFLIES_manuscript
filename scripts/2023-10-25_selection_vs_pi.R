# RRR


# ---- set-up ----

library(data.table)
library(lubridate)

# annotations <- fread(input = "data/2023-10-16_consistent_genes_example_data/combined_file/key_genes_all_genomes-anns_separate_lines.tsv.gz")

genome <- fread("data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")

abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")

tax <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

selection <- fread(file = "data/2023-10-26_selection_per_sample_summaries/all_genome_selection_summaries.tsv.gz")

distances <- fread(file = "data/2023-11-01_multidimensional_SNV_analysis/stats_genomes_combined/all_genomes_all_SNVs_per_sample_distance_stats.tsv.gz")

s.cor <- fread(file = "data/2023-04-09_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-SEASONAL.tsv")
a.cor <- fread(file = "data/2023-04-09_genome_abund-div_correlations/correlations_btwn_diversity_and_abundance-ANNUAL.tsv")

# choose cutoffs

breadth.cutoff <- .5

cor.cutoff <- .3 

coverage.cutoff <- 10

# ---- format files as needed ----

genome[ ,above.breadth.cutoff := (breadth / breadth_expected) >= breadth.cutoff, ]
genome[ ,date := parse_date_time(x = date, orders = "ymd")]

colnames(abund) <- c("sample","genome","abund.perc","Mean.coverage","RPKM","TPM")
abund <- abund[genome != "unmapped"]

tax <- as.data.table(tax)
tax <- tax[winner == TRUE, c("bin.full.name","completeness","contamination","drep.cluster.0.96ANI","length","N50","centrality","score","winner","num.in.cluster",
               "domain","phylum","class","order","family","genus","species")]
colnames(tax)[1] <- "genome"

selection <- selection[ ,c("genome","sample","Npos","Nneg")]
colnames(selection) <- c("genome","sample","N.pos.selected","N.neg.selected")

distances <- distances[ ,-c("dist.obj.order","date","year","yday","season","invasion")]
str(distances) # has column names from how I concatenated it!
distances[ ,`:=`(mean.dist.to.other.samples = as.numeric(mean.dist.to.other.samples),
                 sd.dist.to.other.samples = as.numeric(sd.dist.to.other.samples),
                 mean.dist.to.other.samples.outlier = as.logical(mean.dist.to.other.samples.outlier),
                 mean.dist.to.other.samples.Q4 = as.logical(mean.dist.to.other.samples.Q4),
                 dist.to.dataset.centroid = as.numeric(dist.to.dataset.centroid),
                 dist.to.dataset.centroid.outlier = as.logical(dist.to.dataset.centroid.outlier),
                 dist.to.dataset.centroid.Q4 = as.logical(dist.to.dataset.centroid.Q4))]

s.cor[cor >= cor.cutoff, match := "match"]
s.cor[cor <= -cor.cutoff, match := "offset"]
s.cor[cor > -cor.cutoff & cor < cor.cutoff, match := "other"]
colnames(s.cor)[-1] <- paste0("s.cor.", colnames(s.cor)[-1])

a.cor[cor >= cor.cutoff, match := "match"]
a.cor[cor <= -cor.cutoff, match := "offset"]
a.cor[cor > -cor.cutoff & cor < cor.cutoff, match := "other"]
colnames(a.cor)[-1] <- paste0("a.cor.", colnames(a.cor)[-1])

# ---- merge files ----

combo <- merge(x = abund, y = genome, by = c("genome","sample"), all = T)
combo <- merge(x = combo, y = selection, by = c("genome","sample"), all = T)
combo <- merge(x = combo, y = distances, by = c("genome","sample"), all = T)
combo <- merge(x = combo, y = tax, by = "genome", all = T)
combo <- merge(x = combo, y = s.cor, by = "genome", all = T)
combo <- merge(x = combo, y = a.cor, by = "genome", all = T)

# ---- remove duplicate dates ----

div <- dcast(data = combo, formula = genome ~ sample, value.var = "nucl_diversity")
div <- as.matrix(x = div, rownames = TRUE)

sample.dates <- colnames(div) |>
  substr(start = 3, stop = 12) |>
  parse_date_time(orders = "ymd")

dup.dates <- sample.dates[which(duplicated(sample.dates))]
dup.indexes <- which(sample.dates %in% dup.dates)
dup.samples <- colnames(div)[dup.indexes]

# manually choose which to remove (the prefiltered one, the less standard depth, and all the generous donor samples (there is an actiual sample from same date)):
remove.these.samples <- c("ME2002-07-17pf_3300034102","ME2006-08-18D6_3300036398",
                          "ME2018-11-08GD_3300034116","ME2018-11-08GD_3300042399","ME2018-11-08GD_3300042510","ME2018-11-08GD_3300042901","ME2018-11-08GD_3300042940","ME2018-11-08GD_3300046832")

for (bad.sample in remove.these.samples){
  combo <- combo[sample != bad.sample]
}

# ---- calculate change in nucleotide diversity ----
div <- dcast(data = combo, formula = genome ~ sample, value.var = "nucl_diversity")
div <- as.matrix(x = div, rownames = TRUE)

# check dates are in order
sample.dates <- colnames(div) |>
  substr(start = 3, stop = 12) |>
  parse_date_time(orders = "ymd")
all.equal(order(sample.dates), 1:length(sample.dates)) # TRUE

delta.div <- div[ ,-1] - div[ ,-ncol(div)]

delta.div <- as.data.table(delta.div, keep.rownames = "genome")
delta.div <- melt(data = delta.div, id.vars = "genome", variable.name = "sample", value.name = "delta.nucl_div")

combo <- merge(x = combo, y = delta.div, by = c("genome","sample"), all = TRUE)

# ---- only include samples where genome is present and covered enough to call SNVs ----

combo <- combo[abund.perc > 0 & !is.na(above.breadth.cutoff) & above.breadth.cutoff == TRUE & coverage_median > coverage.cutoff]

# normalize by the organism's max selection & diversity  & SNV distances ----
genome.max.abunds <- combo[ ,.(genome.max.abund = max(abund.perc)), by = .(genome)]
combo <- merge(x = combo, y = genome.max.abunds, by = "genome")
combo[ ,abund.perc.max := abund.perc / genome.max.abund * 100]

genome.max.div <- combo[ ,.(genome.max.nucl_div = max(nucl_diversity, na.rm = T)), by = .(genome)]
combo <- merge(x = combo, y = genome.max.div, by = "genome")
combo[ ,nucl_diversity.perc.max := nucl_diversity / genome.max.nucl_div * 100]

genome.max.npos <- combo[ ,.(genome.max.N.pos.selected = max(N.pos.selected, na.rm = T)), by = .(genome)]
combo <- merge(x = combo, y = genome.max.npos, by = "genome")
combo[ ,N.pos.selected.perc.max := N.pos.selected / genome.max.N.pos.selected * 100]

genome.max.nneg <- combo[ ,.(genome.max.N.neg.selected = max(N.neg.selected, na.rm = T)), by = .(genome)]
combo <- merge(x = combo, y = genome.max.nneg, by = "genome")
combo[ ,N.neg.selected.perc.max := N.neg.selected / genome.max.N.neg.selected * 100]

genome.max.dist <- combo[ ,.(genome.max.dist = max(mean.dist.to.other.samples, na.rm = T)), by = .(genome)]
genome.max.dist[is.infinite(genome.max.dist), genome.max.dist := NA]
combo <- merge(x = combo, y = genome.max.dist, by = "genome")
combo[ ,mean.dist.to.other.samples.perc := mean.dist.to.other.samples / genome.max.dist * 100]


# ---- explore selection vs. nuc.div. ----

# overall, does selection correlate with strain diversity? nope
plot(x = combo$nucl_diversity, y = combo$N.pos.selected, col = adjustcolor("black",alpha.f = .1))

# what about just in more abundant samples?
plot(x = combo[abund.perc >= .01, nucl_diversity], y = combo[abund.perc >= .01, N.pos.selected], col = adjustcolor("black",alpha.f = .1)) # nope
plot(x = combo[abund.perc >= .05, nucl_diversity], y = combo[abund.perc >= .05, N.pos.selected], col = adjustcolor("black",alpha.f = .1)) # nope
# no change now that I added the coverage cutoff, that is a harsher filter.

# what about in only the ones with nuc div patterns?
plot(x = combo[a.cor.match == "match" | a.cor.match == "offset" | a.cor.match == "other" | s.cor.match == "match" | s.cor.match == "offset" | s.cor.match == "other", nucl_diversity], 
     y = combo[a.cor.match == "match" | a.cor.match == "offset" | a.cor.match == "other" | s.cor.match == "match" | s.cor.match == "offset" | s.cor.match == "other", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .1)) # still nope

plot(x = combo[a.cor.match == "offset", nucl_diversity], 
     y = combo[a.cor.match == "offset", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope, but maybe closer

plot(x = combo[s.cor.match == "offset", nucl_diversity], 
     y = combo[s.cor.match == "offset", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[a.cor.match == "match", nucl_diversity], 
     y = combo[a.cor.match == "match", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope, but maybe closer

plot(x = combo[s.cor.match == "match", nucl_diversity], 
     y = combo[s.cor.match == "match", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope, but maybe closer?


# new plots with normalized data ----

# overall, does selection correlate with strain diversity? 
plot(x = combo$nucl_diversity.perc.max, y = combo$N.pos.selected.perc.max, col = adjustcolor("black",alpha.f = .1)) # nope

# what about just in more abundant samples?
plot(x = combo[abund.perc >= .01, nucl_diversity.perc.max], y = combo[abund.perc >= .01, N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .1)) # nope
plot(x = combo[abund.perc >= .05, nucl_diversity.perc.max], y = combo[abund.perc >= .05, N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .1)) # nope

# what about just the more abundant samples *of that genome*- like just the bloom samples?
plot(x = combo[abund.perc.max >= 50, nucl_diversity.perc.max], y = combo[abund.perc.max >= 50, N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .1)) # nope
plot(x = combo[abund.perc.max >= 25, nucl_diversity.perc.max], y = combo[abund.perc.max >= 25, N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .1)) # nope
plot(x = combo[abund.perc.max >= 75, nucl_diversity.perc.max], y = combo[abund.perc.max >= 75, N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .1)) # nope

# what about in only the ones with nuc div patterns?
plot(x = combo[a.cor.match == "match" | a.cor.match == "offset" | a.cor.match == "other" | s.cor.match == "match" | s.cor.match == "offset" | s.cor.match == "other", nucl_diversity.perc.max], 
     y = combo[a.cor.match == "match" | a.cor.match == "offset" | a.cor.match == "other" | s.cor.match == "match" | s.cor.match == "offset" | s.cor.match == "other", N.pos.selected.perc.max], 
     col = adjustcolor("black",alpha.f = .1)) # still nope, but maybe closer?

plot(x = combo[a.cor.match == "offset", nucl_diversity.perc.max], 
     y = combo[a.cor.match == "offset", N.pos.selected.perc.max], 
     col = adjustcolor("black",alpha.f = .4)) # nope

plot(x = combo[s.cor.match == "offset", nucl_diversity.perc.max], 
     y = combo[s.cor.match == "offset", N.pos.selected.perc.max], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[a.cor.match == "match", nucl_diversity.perc.max], 
     y = combo[a.cor.match == "match", N.pos.selected.perc.max], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[s.cor.match == "match", nucl_diversity.perc.max], 
     y = combo[s.cor.match == "match", N.pos.selected.perc.max], 
     col = adjustcolor("black",alpha.f = .4)) # still nope, maybe closer though!

# what about with normalized nuc diversity and not normalized selection?
plot(x = combo[a.cor.match == "offset", nucl_diversity.perc.max], 
     y = combo[a.cor.match == "offset", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # nope

plot(x = combo[s.cor.match == "offset", nucl_diversity.perc.max], 
     y = combo[s.cor.match == "offset", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[a.cor.match == "match", nucl_diversity.perc.max], 
     y = combo[a.cor.match == "match", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[s.cor.match == "match", nucl_diversity.perc.max], 
     y = combo[s.cor.match == "match", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope, maybe closer though!

# what about subsetting to bloom days AND to match ones?
plot(x = combo[a.cor.match == "match" & abund.perc.max >= 50, nucl_diversity.perc.max], 
     y = combo[a.cor.match == "match" & abund.perc.max >= 50, N.pos.selected.perc.max], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[s.cor.match == "match" & abund.perc.max >= 50, nucl_diversity.perc.max], 
     y = combo[s.cor.match == "match" & abund.perc.max >= 50, N.pos.selected.perc.max], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

# what about in the genomes that I THINK this is true in? the ones I've been looking at? ----

acI.B <- combo[genome == "ME2011-09-21_3300043464_group3_bin69"]
acI.C <- combo[genome == "ME2016-07-20_3300033996_group7_bin32"]
acI.A <- combo[genome == "ME2011-09-04_3300044729_group3_bin142"]

plot(x = acI.B[ ,nucl_diversity], y = acI.B[ ,N.pos.selected], col = adjustcolor("black",alpha.f = .5)) # hmmm it's not a linear line
plot(x = acI.B[abund.perc.max >= 25 ,nucl_diversity], y = acI.B[abund.perc.max >= 25 ,N.pos.selected], col = adjustcolor("black",alpha.f = .5)) # it's more like a cutoff. like it has to have high nucl div to see selection but you don't ALWAYS see selection given high nucl diversity

plot(x = acI.C[ ,nucl_diversity], y = acI.C[ ,N.pos.selected], col = adjustcolor("black",alpha.f = .5)) # hmmm it's not a linear line
plot(x = acI.C[abund.perc.max >= 25 ,nucl_diversity], y = acI.C[abund.perc.max >= 25 ,N.pos.selected], col = adjustcolor("black",alpha.f = .5)) # it's more like a cutoff. like it has to have high nucl div to see selection but you don't ALWAYS see selection given high nucl diversity

plot(x = acI.A[ ,nucl_diversity], y = acI.A[ ,N.pos.selected], col = adjustcolor("black",alpha.f = .5)) # but don't see same "cutoff" here, more opposite
plot(x = acI.A[abund.perc.max >= 25 ,nucl_diversity], y = acI.A[abund.perc.max >= 25 ,N.pos.selected], col = adjustcolor("black",alpha.f = .5)) # and no trend at all for the blooms

# If I'm not seeing it this way in the genomes I EXPECT to see it in from looking at them a bunch, 
# what am I doing wrong?

# Maybe it's actually that an INCREASE IN diversity causes selection, but not stable high diversity. 
# Because AS the diversity is increasing, that's when strains are trading off and genes are sweeping.

# Try delta div instead of total div ----

# overall, does selection correlate with strain diversity? 
plot(x = combo$delta.nucl_div, y = combo$N.pos.selected.perc.max, col = adjustcolor("black",alpha.f = .1)) # nope
plot(x = combo$delta.nucl_div, y = combo$N.pos.selected, col = adjustcolor("black",alpha.f = .1)) # nope

# what about just in more abundant samples?
plot(x = combo[abund.perc >= .01, delta.nucl_div], y = combo[abund.perc >= .01, N.pos.selected], col = adjustcolor("black",alpha.f = .1)) # nope
plot(x = combo[abund.perc >= .05, delta.nucl_div], y = combo[abund.perc >= .05, N.pos.selected], col = adjustcolor("black",alpha.f = .1)) # nope

# what about just the more abundant samples *of that genome*- like just the bloom samples?
plot(x = combo[abund.perc.max >= 50, delta.nucl_div], y = combo[abund.perc.max >= 50, N.pos.selected], col = adjustcolor("black",alpha.f = .1)) # nope
plot(x = combo[abund.perc.max >= 25, delta.nucl_div], y = combo[abund.perc.max >= 25, N.pos.selected], col = adjustcolor("black",alpha.f = .1)) # nope
plot(x = combo[abund.perc.max >= 75, delta.nucl_div], y = combo[abund.perc.max >= 75, N.pos.selected], col = adjustcolor("black",alpha.f = .1)) # nope

# what about in only the ones with nuc div patterns?
plot(x = combo[a.cor.match == "match" | a.cor.match == "offset" | a.cor.match == "other" | s.cor.match == "match" | s.cor.match == "offset" | s.cor.match == "other", delta.nucl_div], 
     y = combo[a.cor.match == "match" | a.cor.match == "offset" | a.cor.match == "other" | s.cor.match == "match" | s.cor.match == "offset" | s.cor.match == "other", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .1)) # still nope, but maybe closer?

plot(x = combo[a.cor.match == "offset", delta.nucl_div], 
     y = combo[a.cor.match == "offset", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # nope

plot(x = combo[s.cor.match == "offset", delta.nucl_div], 
     y = combo[s.cor.match == "offset", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[a.cor.match == "match", delta.nucl_div], 
     y = combo[a.cor.match == "match", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[s.cor.match == "match", delta.nucl_div], 
     y = combo[s.cor.match == "match", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope, maybe closer though!

# what about subsetting to bloom days AND to match ones?
plot(x = combo[a.cor.match == "match" & abund.perc.max >= 50, delta.nucl_div], 
     y = combo[a.cor.match == "match" & abund.perc.max >= 50, N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[s.cor.match == "match" & abund.perc.max >= 50, delta.nucl_div], 
     y = combo[s.cor.match == "match" & abund.perc.max >= 50, N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

# what about saying - only within the samples where diversity is changing one way or another
plot(x = combo[s.cor.match == "match" & abund.perc.max >= 50 & delta.nucl_div > 0, delta.nucl_div], 
     y = combo[s.cor.match == "match" & abund.perc.max >= 50 & delta.nucl_div > 0, N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) #hmmm, maybe on the right track?

# do we see it in our genomes this way?
plot(x = acI.B[abund.perc.max >= 25 & delta.nucl_div > 0, delta.nucl_div], y = acI.B[abund.perc.max >= 25 & delta.nucl_div > 0, N.pos.selected], col = adjustcolor("black",alpha.f = .5)) 
plot(x = acI.B[abund.perc.max >= 25 & delta.nucl_div < 0, delta.nucl_div], y = acI.B[abund.perc.max >= 25 & delta.nucl_div < 0, N.pos.selected], col = adjustcolor("black",alpha.f = .5)) 

plot(x = acI.C[abund.perc.max >= 25 & delta.nucl_div > 0, delta.nucl_div], y = acI.C[abund.perc.max >= 25 & delta.nucl_div > 0, N.pos.selected], col = adjustcolor("black",alpha.f = .5)) 
plot(x = acI.C[abund.perc.max >= 25 & delta.nucl_div < 0, delta.nucl_div], y = acI.C[abund.perc.max >= 25 & delta.nucl_div < 0, N.pos.selected], col = adjustcolor("black",alpha.f = .5)) 

plot(x = acI.A[abund.perc.max >= 25 & delta.nucl_div > 0, delta.nucl_div], y = acI.A[abund.perc.max >= 25 & delta.nucl_div > 0, N.pos.selected], col = adjustcolor("black",alpha.f = .5)) 
plot(x = acI.A[abund.perc.max >= 25 & delta.nucl_div < 0, delta.nucl_div], y = acI.A[abund.perc.max >= 25 & delta.nucl_div < 0, N.pos.selected], col = adjustcolor("black",alpha.f = .5))

# color by year
acI.B$color.year <- "green3"
acI.B[year > 2011, color.year := "red2"]
acI.B[year > 2012, color.year := "orange"]
acI.B[year > 2013, color.year := "lightblue4"]

plot(x = acI.B[ ,delta.nucl_div], y = acI.B[ ,N.pos.selected], col = acI.B[ ,color.year]) 
plot(x = acI.B[ ,nucl_diversity], y = acI.B[ ,N.pos.selected], col = acI.B[ ,color.year]) 

color.key <- data.table("season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"),
                        "color" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))
acI.C <- merge(x = acI.C, y = color.key, by = "season")

plot(x = acI.C[ ,delta.nucl_div], y = acI.C[ ,N.pos.selected], col = acI.C[ ,color]) 
plot(x = acI.C[ ,nucl_diversity], y = acI.C[ ,N.pos.selected], col = acI.C[ ,color]) 

# OK. So I see early years, 2012/13, and later years clusters on the acI-B plot
# so the plot looks as expected, but you can't fit a line to it
# I think the conclusion is, it's not as simple as saying more strains == more selection!
# Maybe the nuance is more UNUSUAL strains == more selection, because 2012 and 2013 stood out in SNV frequency space not in nucleotide diversity


# ---- explore selection vs. SNV distances ----

plot(x = combo[ ,mean.dist.to.other.samples], y = combo[ , N.pos.selected], col = adjustcolor("black",alpha.f = .1)) 

plot(x = combo[ ,dist.to.dataset.centroid], y = combo[ , N.pos.selected], col = adjustcolor("black",alpha.f = .1)) 

# these ones do not look great, but at least the acI-B does now
acI.B <- combo[genome == "ME2011-09-21_3300043464_group3_bin69"]
acI.C <- combo[genome == "ME2016-07-20_3300033996_group7_bin32"]
acI.A <- combo[genome == "ME2011-09-04_3300044729_group3_bin142"]

plot(x = acI.B[ ,mean.dist.to.other.samples], y = acI.B[ , N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .8)) 
plot(x = acI.B[ ,dist.to.dataset.centroid], y = acI.B[ , N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .8)) 

plot(x = acI.C[ ,mean.dist.to.other.samples], y = acI.C[ , N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .8)) 
plot(x = acI.C[ ,dist.to.dataset.centroid], y = acI.C[ , N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .8)) 

plot(x = acI.A[ ,mean.dist.to.other.samples], y = acI.A[ , N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .8)) 
plot(x = acI.A[ ,dist.to.dataset.centroid], y = acI.A[ , N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .8)) 

# as percent max?

plot(x = combo[ ,mean.dist.to.other.samples.perc], y = combo[ , N.pos.selected.perc.max], col = adjustcolor("black",alpha.f = .1)) # hard NO
plot(x = combo[ ,mean.dist.to.other.samples], y = combo[ , N.pos.selected], col = adjustcolor("black",alpha.f = .1)) 


# what about in only the ones with nuc div patterns?
plot(x = combo[a.cor.match == "match" | a.cor.match == "offset" | a.cor.match == "other" | s.cor.match == "match" | s.cor.match == "offset" | s.cor.match == "other", mean.dist.to.other.samples], 
     y = combo[a.cor.match == "match" | a.cor.match == "offset" | a.cor.match == "other" | s.cor.match == "match" | s.cor.match == "offset" | s.cor.match == "other", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .1)) # still nope, but maybe closer?

plot(x = combo[a.cor.match == "offset", mean.dist.to.other.samples], 
     y = combo[a.cor.match == "offset", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # nope

plot(x = combo[s.cor.match == "offset", mean.dist.to.other.samples], 
     y = combo[s.cor.match == "offset", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[a.cor.match == "match", mean.dist.to.other.samples], 
     y = combo[a.cor.match == "match", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope

plot(x = combo[s.cor.match == "match", mean.dist.to.other.samples], 
     y = combo[s.cor.match == "match", N.pos.selected], 
     col = adjustcolor("black",alpha.f = .4)) # still nope, maybe closer though!

plot(combo[ ,N.pos.selected] ~ combo[ , mean.dist.to.other.samples], col = adjustcolor("black",alpha.f = .1)) 
x <- lm(combo[ ,N.pos.selected] ~ combo[ , mean.dist.to.other.samples]) 
summary(x)
abline(x)

plot(combo[a.cor.match == "offset" ,N.pos.selected] ~ combo[a.cor.match == "offset" , mean.dist.to.other.samples], col = adjustcolor("black",alpha.f = .5)) 
x <- lm(combo[a.cor.match == "offset" ,N.pos.selected] ~ combo[a.cor.match == "offset" , mean.dist.to.other.samples]) 
summary(x)
abline(x)

plot(combo[s.cor.match == "offset" ,N.pos.selected] ~ combo[s.cor.match == "offset" , mean.dist.to.other.samples], col = adjustcolor("black",alpha.f = .5)) 
x <- lm(combo[s.cor.match == "offset" ,N.pos.selected] ~ combo[s.cor.match == "offset" , mean.dist.to.other.samples]) 
summary(x)
abline(x)

plot(combo[a.cor.match == "match" ,N.pos.selected] ~ combo[a.cor.match == "match" , mean.dist.to.other.samples], col = adjustcolor("black",alpha.f = .5)) 
x <- lm(combo[a.cor.match == "match" ,N.pos.selected] ~ combo[a.cor.match == "match" , mean.dist.to.other.samples]) 
summary(x)
abline(x)

plot(combo[s.cor.match == "match" ,N.pos.selected] ~ combo[s.cor.match == "match" , mean.dist.to.other.samples], col = adjustcolor("black",alpha.f = .5)) 
x <- lm(combo[s.cor.match == "match" ,N.pos.selected] ~ combo[s.cor.match == "match" , mean.dist.to.other.samples]) 
summary(x)
abline(x)
