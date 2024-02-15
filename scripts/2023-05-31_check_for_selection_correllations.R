# RRR

# does amt of selection relate to growth?
# to diversity?


library(data.table)
library(lubridate)

# ---- import ----

# get tax
tax <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
tax <- as.data.table(tax)
tax <- tax[winner == TRUE, .(bin.full.name, domain, phylum, class, order, family, genus, species)]

# get COGs
folder.path <- "data/2023-05-31_COG_summary_data/Selection_Summaries/"
my.files <- list.files(folder.path)

all.genomes.list <- list()
for (f in my.files){
  one.g <- fread(file = file.path(folder.path,f))
  if (nrow(one.g > 0)){ # some had no selection so are empty tables
    one.g$genome <- sub(pattern = "_selected_COGs\\.tsv\\.gz", replacement = "", x = f)
    all.genomes.list <- c(all.genomes.list, list(one.g))
  }
}

all.genomes.table <- rbindlist(l = all.genomes.list)
# rm(all.genomes.list)
all.genomes.table <- merge(x = all.genomes.table, y = tax, by.x = "genome", by.y = "bin.full.name", all.x = T, all.y = F)

# ---- get abund SLOPE ----

abund <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz", sep = "\t")
colnames(abund)[3] <- "rel.abund.perc"
abund <- abund[Genome != "unmapped", .(Genome, Sample,rel.abund.perc)]
abund[ ,date := parse_date_time(substr(Sample, start = 3, stop = 12), orders = "ymd")]
save.abund <- copy(abund)
abund <- dcast(data = abund, formula = Genome ~ Sample, value.var = "rel.abund.perc", fun.aggregate = mean)
abund[1:5,1:5]

abund.1 <- as.matrix(x = abund, rownames = T)
abund.2 <- abund.1
abund.1 <- abund.1[ ,-1]
abund.2 <- abund.2[ ,-ncol(abund.2)]
abund.1[1:5, 1:5]
abund.2[1:5,1:5]
abund.1[1:5, (ncol(abund.1)-5):ncol(abund.1)]
abund.2[1:5, (ncol(abund.2)-5):ncol(abund.2)]

change <- abund.1 - abund.2
change[1:5,1:5]

dates.1 <- parse_date_time(x = substr(colnames(abund.1),start = 3, stop = 12), "ymd")
dates.2 <- parse_date_time(x = substr(colnames(abund.2),start = 3, stop = 12), "ymd")
dates.1[1:5]
dates.2[1:5]
days.apart <- as.numeric(dates.1 - dates.2)
dim(change)

x <- matrix(1:12,4,3)
y <- 1:4
z <- 1:3
x/y
x/z
# math is applied down the columns

change <- t(change)
change[1:5,1:5]
growth.rate <- change / days.apart

growth.rate <- as.data.table(growth.rate, keep.rownames = "Sample")
growth.rate <- melt(data = growth.rate, id.vars = "Sample", variable.name = "Genome", value.name = "growth.rate")

growth.rate <- merge(x = growth.rate, y = save.abund, by = c("Genome","Sample"))

# growth rate units: rel abund change per day

all.genomes.table <- merge(x = all.genomes.table[ ,-c("date")], y = growth.rate[ ,-c("date")], by.x = c("sample","genome"), by.y = c("Sample","Genome"), all = F)

# plot(rel.abund.perc ~ Npos, data = all.genomes.table, col = adjustcolor("red3",.1)) # anti-correlated actually!
# plot(rel.abund.perc ~ Nneg, data = all.genomes.table, col = adjustcolor("blue2",.1))
# plot(rel.abund.perc ~ Npos.perc, data = all.genomes.table, col = adjustcolor("red3",.1)) # anti-correlated actually!
# plot(rel.abund.perc ~ Nneg, data = all.genomes.table, col = adjustcolor("blue2",.1))
# 
# plot(growth.rate ~ Npos, data = all.genomes.table, col = adjustcolor("red3",.1)) # anti-correlated actually!
# plot(growth.rate ~ Nneg, data = all.genomes.table, col = adjustcolor("blue2",.1))
# plot(growth.rate ~ Npos.perc, data = all.genomes.table, col = adjustcolor("red3",.1)) # anti-correlated actually!
# plot(growth.rate ~ Nneg, data = all.genomes.table, col = adjustcolor("blue2",.1))
# 
# plot(growth.rate ~ Npos.perc, data = all.genomes.table[genome == unique(genome)[300]], col = adjustcolor("red3",.5)) # anti-corr
# plot(growth.rate ~ Nneg, data = all.genomes.table[genome == unique(genome)[300]], col = adjustcolor("blue2",.1))
# 
# plot(abs(growth.rate) ~ Npos.perc, data = all.genomes.table, col = adjustcolor("red3",.5)) # anti-corr
# plot(abs(growth.rate) ~ Nneg, data = all.genomes.table, col = adjustcolor("blue2",.1))

# OK so that was a bust...


library(limony)
data("limony")
data("key")
data("seasons")

library(vegan)

me <- limony$av$seqID

pres.abs <- me > 0
richness <- colSums(pres.abs)

shannon <- apply(X = me, MARGIN = 2, FUN = diversity, index = "shannon")

all.equal(names(richness), key$in.R.colnames)
all.equal(names(shannon), key$in.R.colnames)
key <- cbind(key, "Richness" = richness, "Shannon.Diversity" = shannon)
head(key)
which(duplicated(key$Sample.Date)) # pf ones
key <- aggregate(x = key[ ,c("Richness","Shannon.Diversity")], by = list(key$Sample.Date), FUN = mean)

all.genomes.table[ , date := parse_date_time(substr(sample, start = 3, stop = 12), orders = "ymd", tz = "Etc/GMT-5")]
genomes <- all.genomes.table[ ,.(Npos = mean(Npos, na.rm = T), Nneg = mean(Nneg, na.rm = T), Npos.perc = mean(Npos.perc, na.rm = T), Nneg.perc = mean(Nneg.perc, na.rm = T)), by = .(COG_category, sample, growth.rate, rel.abund.perc, date)]
genomes <- all.genomes.table[ ,.(Npos = sum(Npos, na.rm = T), Nneg = sum(Nneg, na.rm = T), Npos.perc = sum(Npos.perc, na.rm = T), Nneg.perc = sum(Nneg.perc, na.rm = T)), by = .(date, growth.rate, rel.abund.perc)]

## same
# plot(rel.abund.perc ~ Npos, data = genomes, col = adjustcolor("red3",.1)) # anti-correlated actually!
# plot(rel.abund.perc ~ Nneg, data = genomes, col = adjustcolor("blue2",.1))
# plot(rel.abund.perc ~ Npos.perc, data = genomes, col = adjustcolor("red3",.1)) # anti-correlated actually!
# plot(rel.abund.perc ~ Nneg, data = genomes, col = adjustcolor("blue2",.1))
# 
# plot(growth.rate ~ Npos, data = genomes, col = adjustcolor("red3",.1)) # anti-correlated actually!
# plot(growth.rate ~ Nneg, data = genomes, col = adjustcolor("blue2",.1))
# plot(growth.rate ~ Npos.perc, data = genomes, col = adjustcolor("red3",.1)) # anti-correlated actually!
# plot(growth.rate ~ Nneg, data = genomes, col = adjustcolor("blue2",.1))
# 
# plot(growth.rate ~ Npos.perc, data = genomes[genome == unique(genome)[300]], col = adjustcolor("red3",.5)) # anti-corr
# plot(growth.rate ~ Nneg, data = genomes[genome == unique(genome)[300]], col = adjustcolor("blue2",.1))
# 
# plot(abs(growth.rate) ~ Npos.perc, data = genomes, col = adjustcolor("red3",.5)) # anti-corr
# plot(abs(growth.rate) ~ Nneg, data = genomes, col = adjustcolor("blue2",.1))


genomes <- merge(x = genomes, y = key, by.x = "date", by.y = "Group.1", all = F)

plot(Richness ~ Npos, data = genomes)
plot(Shannon.Diversity ~ Npos, data = genomes)
plot(Richness ~ Nneg, data = genomes)
plot(Shannon.Diversity ~ Nneg, data = genomes)
