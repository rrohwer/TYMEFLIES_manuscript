# RRR
# In which samples are the SNVs outliers?
#   --> calc distance from each sample to all other samples
#   --> calc distance from each sample to the centroid of all samples
#   --> identify outlier and Q4 samples
# Are years different in their SNV compositions? 
#   --> calc distance of each year centroid from all other year centroids
#   --> permanova to see if year centroids are different from each other
#   --> permdisp to see if year dispersions are different from each other
# Is there long-term drift?
#   --> calc distance of each year from the first year
# Are seasons different in their SNV compositions? 
#   --> calc distance of each season centroid from all other season centroids
#   --> permanova to see if season centroids are different from each other
#   --> permdisp to see if season dispersions are different from each other
#   --> IF seasons are different, revisit the "reset" idea
# How widespread is having 2012/2013 as an outlier? (subset of 2012, just the samples after the jump)
#   --> calc distance of 2012/2013 centroid to pre, post, and dataset centroids
#   --> permanova and permdisp to see if it's significantly different from them
#   --> compare distance of 2012/2013 to all other years vs those years to each other

# note: working with dist objects, it is not a simple matrix
# https://cran.r-project.org/web/packages/usedist/readme/README.html#:~:text=In%20R%2C%20the%20dist(),it%20look%20like%20a%20matrix.
# wow the usedist package is simple and seems to do everything I want besides the stats calcs!

# note: adonis2 is vegan's implementation of a PERMANOVA test. syntax:
# permanova.years <- adonis2(dist.obj ~ year, data = sample.key, permutations = 999)
# permanova.years$`Pr(>F)`[1] # the p-value
# but this gives a single p-value, so need to make it do pairwise comparisons yourself

# note: betadisper is vegan's implementation of a PERMDISP test. syntax:
# permdisp.years <- betadisper(d = dist.obj, group = sample.key$year)
# x <- permutest(permdisp.years) # this gives you a single p-value


# ---- set up ----

library(data.table)
library(lubridate)
library(vegan)
library(usedist)

userprefs <- commandArgs(trailingOnly = TRUE)
dist.file <- userprefs[1]
sample.key <- userprefs[2]
output.per.sample.stats <- userprefs[3]
output.per.year.stats <- userprefs[4]
output.pairwise.year.stats <- userprefs[5]
output.per.season.stats <- userprefs[6]
output.pairwise.season.stats <- userprefs[7]
threads <- as.numeric(userprefs[8])

# # local paths testing
# dist.file <- "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-21_3300043464_group3_bin69_all_SNV_euclidean_distance_matrix.rds"
# sample.key <- "data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv"
# output.per.sample.stats <- "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-21_3300043464_group3_bin69_all_SNV_per_sample_SNV_stats.tsv.gz"
# output.per.year.stats <- "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-21_3300043464_group3_bin69_all_SNV_per_year_SNV_stats.tsv.gz"
# output.pairwise.year.stats <- "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-21_3300043464_group3_bin69_all_SNV_pairwise_year_SNV_stats.tsv.gz"
# output.per.season.stats <- "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-21_3300043464_group3_bin69_all_SNV_per_season_SNV_stats.tsv.gz"
# output.pairwise.season.stats <- "data/2023-11-01_multidimensional_SNV_analysis/ME2011-09-21_3300043464_group3_bin69_all_SNV_pairwise_season_SNV_stats.tsv.gz"
# threads <- 1


# ---- import and format ----

# not all genomes HAVE a dist object, they had to have enough coverage. skip if missing the file:
if (!file.exists(dist.file)){
  quit(save = "no", status = 0)
}

my.genome <- sub(pattern = "^.*ME", replacement = "", x = dist.file)
my.genome <- sub(pattern = "_all_SNV_euclidean_distance_matrix.rds", replacement = "", x = my.genome)
my.genome <- sub(pattern = "_nonsynonymous_SNV_euclidean_distance_matrix.rds", replacement = "", x = my.genome)
my.genome <- paste0("ME",my.genome)
cat("Processing", my.genome,"\n")

dist.obj <- readRDS(dist.file)

sample.key <- fread(file = sample.key, colClasses = "character")

# make sample.key match the order of samples in dist.obj
dist.order <- data.table("sample" = attributes(dist.obj)$Labels)
dist.order[ ,dist.obj.order := 1:nrow(dist.order)]
sample.key <- merge(x = dist.order, y = sample.key, by = "sample", all.x = TRUE, all.y = FALSE)
sample.key <- sample.key[order(dist.obj.order)]


# ---- calc distance from each sample to all other samples ----

dist.mat <- as.matrix(dist.obj)

# make diagonals NA
for (r in 1:nrow(dist.mat)){ 
  dist.mat[r,r] <- NA
}

# calculate stats
dist.dt <- as.data.table(dist.mat, keep.rownames = "sample")
dist.long <- melt(data = dist.dt, id.vars = "sample", value.name = "dist", variable.name = "sample.2")
dist.long <- dist.long[!is.na(dist)]

dist.stats <- dist.long[ ,.(mean.dist.to.other.samples = mean(dist),
                            sd.dist.to.other.samples = sd(dist)), by = sample]

# anova tells you if any group is different
# Tukey HSD test tells you which groups are different
# my.anova <- aov(dist ~ sample, data = dist.long)
# TukeyHSD(my.anova)
# But this is pairwise, sample-to-sample
# Could tally up how many samples it's different from

# seems simpler to just say which mean distances are outliers?
box.stats <- boxplot(dist.stats[ ,mean.dist.to.other.samples], plot = F)
outlier.cutoff <- min(box.stats$out[box.stats$out > box.stats$stats[5,1]]) # use >=
Q.stats <- quantile(dist.stats[ ,mean.dist.to.other.samples])
Q4.cutoff <- Q.stats[4] # use >

dist.stats[ ,`:=`(mean.dist.to.other.samples.outlier = mean.dist.to.other.samples >= outlier.cutoff,
                  mean.dist.to.other.samples.Q4 = mean.dist.to.other.samples > Q4.cutoff)]


# ---- calc distance from each sample to the centroid of all samples ----

data.centroid <- dist_to_centroids(d = dist.obj, g = rep("All",nrow(dist.stats)))
colnames(data.centroid) <- c("sample","group","dist.to.dataset.centroid")
dist.stats <- merge(x = dist.stats, y = data.centroid[ ,c(1,3)], by = "sample")

# which centroid distances are outliers?
box.stats <- boxplot(dist.stats[ ,dist.to.dataset.centroid], plot = F)
outlier.cutoff <- min(box.stats$out[box.stats$out > box.stats$stats[5,1]]) # use >=
Q.stats <- quantile(dist.stats[ ,dist.to.dataset.centroid])
Q4.cutoff <- Q.stats[4] # use >

dist.stats[ ,`:=`(dist.to.dataset.centroid.outlier = dist.to.dataset.centroid >= outlier.cutoff,
                  dist.to.dataset.centroid.Q4 = dist.to.dataset.centroid > Q4.cutoff)]


# ---- YEARS ----

# ---- calc distance of each year centroid from all other year centroids ----

dist.stats <- merge(x = sample.key, y = dist.stats, by = "sample")

# get a list of samples in each year
year.list <- list()
for (y in unique(dist.stats$year)){
  year.samples <- dist.stats[year == y, sample]
  year.samples <- list(year.samples)
  names(year.samples) <- y
  year.list <- c(year.list, year.samples)
}

# make a year-to-year comparison table
dist.years.mat <- matrix(data = as.numeric(NA), nrow = length(unique(dist.stats$year)), ncol = length(unique(dist.stats$year)))
row.names(dist.years.mat) <- unique(dist.stats$year)
colnames(dist.years.mat) <- row.names(dist.years.mat)
dist.years.dt <- as.data.table(dist.years.mat, keep.rownames = "year.1")
dist.years.dt <- melt(data = dist.years.dt, id.vars = "year.1", variable.name = "year.2", value.name = "dist")
dist.years.dt[ ,year.2 := as.character(year.2)]
dist.years.dt <- dist.years.dt[year.1 != year.2] # remove "diagonal", so skip comparing ie 2000v2000
# don't do this next line. it would be more efficient, but then collapsing to per year stats gets to confusing so just let it take a little longer
# dist.years.dt <- dist.years.dt[year.1 > year.2] # remove "upper" part of dist matrix, i.e. keep only unique pairs, so 2000v2001 == 2001v2000

# calculate distances btwn year centroids
for (r in 1:nrow(dist.years.dt)){
  dist.years.dt[r,dist := dist_between_centroids(d = dist.obj, idx1 = year.list[[dist.years.dt[r,year.1]]], idx2 = year.list[[dist.years.dt[r,year.2]]])]
}

# are these centroids significantly different?
for (r in 1:nrow(dist.years.dt)){
  # subset distance matrix and sample key to the two years in question, keep name order matching
  my.key <- sample.key[year == dist.years.dt[r, year.1] | year == dist.years.dt[r, year.2]]
  my.dist <- dist_subset(d = dist.obj, idx = my.key$sample)
  
  pairwise.permanova <- adonis2(my.dist ~ year, data = my.key, permutations = 999)
  my.pval <- pairwise.permanova$`Pr(>F)`[1]
  
  dist.years.dt[r, permanova.pval := my.pval]
}

# make a by-year stats summary table
dist.years.stats <- dist.years.dt[ , .(mean.dist.to.other.years = mean(dist),
                                       sd.dist.to.other.years = sd(dist)), by = year.1]

box.stats <- boxplot(dist.years.stats[ ,mean.dist.to.other.years], plot = F)
outlier.cutoff <- min(box.stats$out[box.stats$out > box.stats$stats[5,1]]) # use >=
Q.stats <- quantile(dist.years.stats[ ,mean.dist.to.other.years])
Q4.cutoff <- Q.stats[4] # use >

dist.years.stats[ ,`:=`(mean.dist.to.other.years.outlier = mean.dist.to.other.years >= outlier.cutoff,
                  mean.dist.to.other.years.Q4 = mean.dist.to.other.years > Q4.cutoff)]


# ---- calc distance from each year centroid to the centroid of all other samples ----

for(r in 1:nrow(dist.years.stats)){ # don't understand why this loop is necessary but get an error if simply 1 line command with idx1 = year.list[[year.1]]
  dist.years.stats[r,dist.to.dataset.centroid := dist_between_centroids(d = dist.obj, 
                                                                        idx1 = year.list[[dist.years.stats[r,year.1]]], 
                                                                        idx2 = dist.stats[!(sample %in% year.list[[dist.years.stats[r,year.1]]]), sample])]
} 

# which centroid distances are outliers?
box.stats <- boxplot(dist.years.stats[ ,dist.to.dataset.centroid], plot = F)
outlier.cutoff <- min(box.stats$out[box.stats$out > box.stats$stats[5,1]]) # use >=
Q.stats <- quantile(dist.years.stats[ ,dist.to.dataset.centroid])
Q4.cutoff <- Q.stats[4] # use >

dist.years.stats[ ,`:=`(dist.to.dataset.centroid.outlier = dist.to.dataset.centroid >= outlier.cutoff,
                  dist.to.dataset.centroid.Q4 = dist.to.dataset.centroid > Q4.cutoff)]

# are year centroids significantly different from centroids of all other years?
for (r in 1:nrow(dist.years.stats)){
  # add new groups of year, not-year to key
  my.key <- copy(sample.key)
  my.key[year == dist.years.stats[r, year.1], year.group := "the.year"]
  my.key[year != dist.years.stats[r, year.1], year.group := "rest"]
  
  pairwise.permanova <- adonis2(dist.obj ~ year.group, data = my.key, permutations = 999)
  my.pval <- pairwise.permanova$`Pr(>F)`[1]
  
  dist.years.stats[r, dist.to.dataset.centroid.permanova.pval := my.pval]
}


# ---- add distance from each year to the first year ----

directional.change <- dist.years.dt[year.2 == min(year.2)]
colnames(directional.change)[colnames(directional.change) == "dist"] <- "dist.to.1st.year.centroid"
colnames(directional.change)[colnames(directional.change) == "permanova.pval"] <- "dist.to.1st.year.centroid.permanova.pval"
directional.change <- directional.change[ ,-c("year.2")]

dist.years.stats <- merge(x = dist.years.stats, y = directional.change, by = "year.1", all = T)

# is there a slope of increasing distance?
# plot(dist.to.1st.year.centroid ~ year.1, data = dist.years.stats)
# x <- lm(dist.to.1st.year.centroid ~ as.numeric(year.1), data = dist.years.stats)
# y <- summary(x)
# y$adj.r.squared
# this doesn't fit in the table well, incorporate it later


# ---- add dispersion of each year ----

permdisp.years <- betadisper(d = dist.obj, group = sample.key$year)

# add pairwise p-values to dist.years.dt based on Tukey multiple comparison of means
pairwise.permdisp <- TukeyHSD(x = permdisp.years)
pairwise.permdisp <- pairwise.permdisp$group 
pairwise.permdisp <- as.data.table(pairwise.permdisp, keep.rownames = "pair")
pairwise.permdisp[ ,`:=`(year.1 = substr(x = pair, start = 1, stop = 4),
                         year.2 = substr(x = pair, start = 6, stop = 9))]
pairwise.permdisp <- pairwise.permdisp[ ,.(year.1, year.2, "permdisp.tukeyHSD.pval" = `p adj`)]

# get the "upper" part of the distance matrix added in too, b/c the longer but repetitive table will make plotting things and stats easier down the road (ie group by year.1 will include all comparisons)
pairwise.permdisp <- dcast(data = pairwise.permdisp, formula = year.1 ~ year.2, value.var = "permdisp.tukeyHSD.pval")
pairwise.permdisp <- as.matrix(x = pairwise.permdisp, rownames = T)
pairwise.permdisp <- rbind(NA, pairwise.permdisp)
pairwise.permdisp <- cbind(pairwise.permdisp, NA)
row.names(pairwise.permdisp)[1] <- colnames(pairwise.permdisp)[1]
colnames(pairwise.permdisp)[ncol(pairwise.permdisp)] <- row.names(pairwise.permdisp)[nrow(pairwise.permdisp)]
pairwise.permdisp <- as.dist(m = pairwise.permdisp, diag = F, upper = F)
pairwise.permdisp <- as.matrix(pairwise.permdisp)
pairwise.permdisp <- as.data.table(x = pairwise.permdisp, keep.rownames = "year.1")
pairwise.permdisp <- melt(data = pairwise.permdisp, id.vars = "year.1", variable.name = "year.2", value.name = "permdisp.tukeyHSD.pval")
pairwise.permdisp <- pairwise.permdisp[year.1 != year.2]

# now can merge and keep all the duplicate pairings
dist.years.dt <- merge(x = dist.years.dt, y = pairwise.permdisp, by = c("year.1","year.2"))

# add outlier dispersion info to dist.years.stats table
permdisp.years <- data.frame("year.1" = names(permdisp.years$group.distances), "year.centroid.dispersion" = permdisp.years$group.distances) |>
  as.data.table()

box.stats <- boxplot(permdisp.years[ ,year.centroid.dispersion], plot = F)
outlier.cutoff <- min(box.stats$out[box.stats$out > box.stats$stats[5,1]]) # use >=
Q.stats <- quantile(dist.years.stats[ ,dist.to.dataset.centroid])
Q4.cutoff <- Q.stats[4] # use >

permdisp.years[ ,`:=`(year.centroid.dispersion.outlier = year.centroid.dispersion >= outlier.cutoff,
                      year.centroid.dispersion.Q4 = year.centroid.dispersion > Q4.cutoff)]

dist.years.stats <- merge(x = dist.years.stats, y = permdisp.years, by = "year.1")



# ---- SEASONS ----

# ---- calc distance of each season centroid from all other season centroids ----

# get a list of samples in each season
season.list <- list()
for (s in unique(dist.stats$season)){
  season.samples <- dist.stats[season == s, sample]
  season.samples <- list(season.samples)
  names(season.samples) <- s
  season.list <- c(season.list, season.samples)
}

# make a season-to-season comparison table
dist.seasons.mat <- matrix(data = as.numeric(NA), nrow = length(unique(dist.stats$season)), ncol = length(unique(dist.stats$season)))
row.names(dist.seasons.mat) <- unique(dist.stats$season)
colnames(dist.seasons.mat) <- row.names(dist.seasons.mat)
dist.seasons.dt <- as.data.table(dist.seasons.mat, keep.rownames = "season.1")
dist.seasons.dt <- melt(data = dist.seasons.dt, id.vars = "season.1", variable.name = "season.2", value.name = "dist")
dist.seasons.dt[ ,season.2 := as.character(season.2)]
dist.seasons.dt <- dist.seasons.dt[season.1 != season.2] # remove "diagonal", so skip comparing ie 2000v2000
# don't do this next line. it would be more efficient, but then collapsing to per season stats gets to confusing so just let it take a little longer
# dist.seasons.dt <- dist.seasons.dt[season.1 > season.2] # remove "upper" part of dist matrix, i.e. keep only unique pairs, so 2000v2001 == 2001v2000

# calculate distances btwn season centroids
for (r in 1:nrow(dist.seasons.dt)){
  dist.seasons.dt[r,dist := dist_between_centroids(d = dist.obj, idx1 = season.list[[dist.seasons.dt[r,season.1]]], idx2 = season.list[[dist.seasons.dt[r,season.2]]])]
}

# are these centroids significantly different?
for (r in 1:nrow(dist.seasons.dt)){
  # subset distance matrix and sample key to the two seasons in question, keep name order matching
  my.key <- sample.key[season == dist.seasons.dt[r, season.1] | season == dist.seasons.dt[r, season.2]]
  my.dist <- dist_subset(d = dist.obj, idx = my.key$sample)
  
  pairwise.permanova <- adonis2(my.dist ~ season, data = my.key, permutations = 999)
  my.pval <- pairwise.permanova$`Pr(>F)`[1]
  
  dist.seasons.dt[r, permanova.pval := my.pval]
}

# make a by-season stats summary table
dist.seasons.stats <- dist.seasons.dt[ , .(mean.dist.to.other.seasons = mean(dist),
                                       sd.dist.to.other.seasons = sd(dist)), by = season.1]

box.stats <- boxplot(dist.seasons.stats[ ,mean.dist.to.other.seasons], plot = F)
outlier.cutoff <- min(box.stats$out[box.stats$out > box.stats$stats[5,1]]) # use >=
Q.stats <- quantile(dist.seasons.stats[ ,mean.dist.to.other.seasons])
Q4.cutoff <- Q.stats[4] # use >

dist.seasons.stats[ ,`:=`(mean.dist.to.other.seasons.outlier = mean.dist.to.other.seasons >= outlier.cutoff,
                        mean.dist.to.other.seasons.Q4 = mean.dist.to.other.seasons > Q4.cutoff)]


# ---- calc distance from each season centroid to the centroid of all other samples ----

for(r in 1:nrow(dist.seasons.stats)){ # don't understand why this loop is necessary but get an error if simply 1 line command with idx1 = season.list[[season.1]]
  dist.seasons.stats[r,dist.to.dataset.centroid := dist_between_centroids(d = dist.obj, 
                                                                        idx1 = season.list[[dist.seasons.stats[r,season.1]]], 
                                                                        idx2 = dist.stats[!(sample %in% season.list[[dist.seasons.stats[r,season.1]]]), sample])]
} 

# which centroid distances are outliers?
box.stats <- boxplot(dist.seasons.stats[ ,dist.to.dataset.centroid], plot = F)
outlier.cutoff <- min(box.stats$out[box.stats$out > box.stats$stats[5,1]]) # use >=
Q.stats <- quantile(dist.seasons.stats[ ,dist.to.dataset.centroid])
Q4.cutoff <- Q.stats[4] # use >

dist.seasons.stats[ ,`:=`(dist.to.dataset.centroid.outlier = dist.to.dataset.centroid >= outlier.cutoff,
                        dist.to.dataset.centroid.Q4 = dist.to.dataset.centroid > Q4.cutoff)]

# are season centroids significantly different from centroids of all other seasons?
for (r in 1:nrow(dist.seasons.stats)){
  # add new groups of season, not-season to key
  my.key <- copy(sample.key)
  my.key[season == dist.seasons.stats[r, season.1], season.group := "the.season"]
  my.key[season != dist.seasons.stats[r, season.1], season.group := "rest"]
  
  pairwise.permanova <- adonis2(dist.obj ~ season.group, data = my.key, permutations = 999)
  my.pval <- pairwise.permanova$`Pr(>F)`[1]
  
  dist.seasons.stats[r, dist.to.dataset.centroid.permanova.pval := my.pval]
}


# ---- add distance from each season to the first season ----

directional.change <- dist.seasons.dt[season.2 == min(season.2)]
colnames(directional.change)[colnames(directional.change) == "dist"] <- "dist.to.1st.season.centroid"
colnames(directional.change)[colnames(directional.change) == "permanova.pval"] <- "dist.to.1st.season.centroid.permanova.pval"
directional.change <- directional.change[ ,-c("season.2")]

dist.seasons.stats <- merge(x = dist.seasons.stats, y = directional.change, by = "season.1", all = T)


# ---- add dispersion of each season ----

permdisp.seasons <- betadisper(d = dist.obj, group = sample.key$season)

# add pairwise p-values to dist.seasons.dt based on Tukey multiple comparison of means
pairwise.permdisp <- TukeyHSD(x = permdisp.seasons)
pairwise.permdisp <- pairwise.permdisp$group 
pairwise.permdisp <- as.data.table(pairwise.permdisp, keep.rownames = "pair")
pairwise.permdisp[ ,`:=`(season.1 = sub(pattern = "-.*$", replacement = "", x = pair),
                         season.2 = sub(pattern = "^.*-", replacement = "", x = pair))]
pairwise.permdisp <- pairwise.permdisp[ ,.(season.1, season.2, "permdisp.tukeyHSD.pval" = `p adj`)]

# get the "upper" part of the distance matrix added in too, b/c the longer but repetitive table will make plotting things and stats easier down the road (ie group by season.1 will include all comparisons)
pairwise.permdisp <- dcast(data = pairwise.permdisp, formula = season.1 ~ season.2, value.var = "permdisp.tukeyHSD.pval")
pairwise.permdisp <- as.matrix(x = pairwise.permdisp, rownames = T)
pairwise.permdisp <- rbind(NA, pairwise.permdisp)
pairwise.permdisp <- cbind(pairwise.permdisp, NA)
row.names(pairwise.permdisp)[1] <- colnames(pairwise.permdisp)[1]
colnames(pairwise.permdisp)[ncol(pairwise.permdisp)] <- row.names(pairwise.permdisp)[nrow(pairwise.permdisp)]
pairwise.permdisp <- as.dist(m = pairwise.permdisp, diag = F, upper = F)
pairwise.permdisp <- as.matrix(pairwise.permdisp)
pairwise.permdisp <- as.data.table(x = pairwise.permdisp, keep.rownames = "season.1")
pairwise.permdisp <- melt(data = pairwise.permdisp, id.vars = "season.1", variable.name = "season.2", value.name = "permdisp.tukeyHSD.pval")
pairwise.permdisp <- pairwise.permdisp[season.1 != season.2]

# now can merge and keep all the duplicate pairings
dist.seasons.dt <- merge(x = dist.seasons.dt, y = pairwise.permdisp, by = c("season.1","season.2"))

# add outlier dispersion info to dist.seasons.stats table
permdisp.seasons <- data.frame("season.1" = names(permdisp.seasons$group.distances), "season.centroid.dispersion" = permdisp.seasons$group.distances) |>
  as.data.table()

box.stats <- boxplot(permdisp.seasons[ ,season.centroid.dispersion], plot = F)
outlier.cutoff <- min(box.stats$out[box.stats$out > box.stats$stats[5,1]]) # use >=
Q.stats <- quantile(dist.seasons.stats[ ,dist.to.dataset.centroid])
Q4.cutoff <- Q.stats[4] # use >

permdisp.seasons[ ,`:=`(season.centroid.dispersion.outlier = season.centroid.dispersion >= outlier.cutoff,
                      season.centroid.dispersion.Q4 = season.centroid.dispersion > Q4.cutoff)]

dist.seasons.stats <- merge(x = dist.seasons.stats, y = permdisp.seasons, by = "season.1")



# ---- 2012 STEP CHANGE ----

# ---- calc distance of 2012/2013 affected samples from all others

# It's not clear where to draw the cutoff.
# nucl div is up > 6/2/2012
# NMDS is jumped > 8/3/2012


# ---- EXPORT DATA ----


# ---- add genome column so tables can be easily concatenated

dist.stats <- cbind(data.table("genome" = my.genome), dist.stats)
dist.years.dt <- cbind(data.table("genome" = my.genome), dist.years.dt)
dist.years.stats <- cbind(data.table("genome" = my.genome), dist.years.stats)
dist.seasons.dt <- cbind(data.table("genome" = my.genome), dist.seasons.dt)
dist.seasons.stats <- cbind(data.table("genome" = my.genome), dist.seasons.stats)


# ---- export files ---- 

# by-sample info
# dist.long # no need to save, basically just the dist.obj
fwrite(x = dist.stats, file = output.per.sample.stats, sep = "\t", nThread = threads)

# by-year info
fwrite(x = dist.years.dt, file = output.pairwise.year.stats, sep = "\t", nThread = threads)
fwrite(x = dist.years.stats, file = output.per.year.stats, sep = "\t", nThread = threads)

# by-season info
fwrite(x = dist.seasons.dt, file = output.pairwise.season.stats, sep = "\t", nThread = threads)
fwrite(x = dist.seasons.stats, file = output.per.season.stats, sep = "\t", nThread = threads)


# ~end~




