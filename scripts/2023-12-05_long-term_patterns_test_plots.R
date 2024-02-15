

library(data.table)
library(lubridate)

snv.stats <- fread(file = "data/2023-12-04_long-term_change_classifier/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")
tax <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
genome.info <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
abunds <- fread(file = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
sample.key <- fread("data/2023-11-01_multidimensional_SNV_analysis/sample_key.tsv", colClasses = c("date" = "character"))

# ---- combine tables ----

# get number of times a genome occurs, as well as other abundance stats

colnames(abunds)[3] <- "rel.abund"
abunds <- abunds[ ,1:3]

sample.key[ ,date := parse_date_time(date, orders = "ymd")]

abunds <- merge(x = sample.key, y = abunds, by.x = "sample", by.y = "Sample", all = TRUE)

# combine with genome info to include coverage info
abunds <- merge(x = abunds, y = genome.info[ ,.(genome,coverage,breadth,breadth_expected,coverage_median,sample)], by.x = c("sample","Genome"), by.y = c("sample","genome"), all = T)

# remove duplicate samples
# manually choose which to remove (the prefiltered one, the less standard depth, and all the generous donor samples (there is an actual sample from same date)):
remove.these.samples <- c("ME2002-07-17pf_3300034102","ME2006-08-18D6_3300036398",
                          "ME2018-11-08GD_3300034116","ME2018-11-08GD_3300042399","ME2018-11-08GD_3300042510","ME2018-11-08GD_3300042901","ME2018-11-08GD_3300042940","ME2018-11-08GD_3300046832")
sample.key <- sample.key[!(sample %in% remove.these.samples)]
any(duplicated(sample.key$date))

abunds <- merge(x = sample.key, y = abunds, by = "sample", all.x = TRUE, all.y = FALSE)

abunds[ ,above.breadth.cutoff := (breadth / breadth_expected) >= .5]
abunds$present <- FALSE
abunds[above.breadth.cutoff == TRUE & rel.abund > 0 ,present := TRUE]
abunds$present.med.cov.10 <- FALSE
abunds[above.breadth.cutoff == TRUE & rel.abund > 0 & coverage_median > 10, present.med.cov.10 := TRUE]

abunds <- abunds[ , .(times.present = sum(present),
                      times.present.med.cov.10 = sum(present.med.cov.10),
                      sum.abund = sum(rel.abund),
                      mean.abund = mean(rel.abund),
                      max.abund = max(rel.abund),
                      med.abund = median(rel.abund)), by = .(Genome)]

tax <- merge(x = abunds, y = tax, by.x = "Genome", by.y = "bin.full.name")

snv.stats <- merge(x = snv.stats, y = tax, by.x = "genome", by.y = "Genome")

snv.stats[is.na(phylum), phylum := "unknown"]

snv.stats[ ,`:=`(max.disturbance.start.date = parse_date_time(max.disturbance.start.date, orders = "ymd"),
                 breakpoint.date = parse_date_time(breakpoint.date, orders = "ymd"))]

summary(snv.stats$times.present)
summary(snv.stats$times.present.med.cov.10)
summary(snv.stats$mean.abund)
nrow(snv.stats)

# ---- classify as seasonal or not ---- 

nrow(snv.stats[peak.1.is.year == T | peak.2.is.year == T | peak.3.is.year == T | peak.4.is.year == T | peak.5.is.year])
# oops need to fix this

# ---- pull out the year of abrupt changes ----

snv.changes <- snv.stats[Classified.LT.Change == "step" | Classified.LT.Change == "disturbance"]
snv.changes[Classified.LT.Change == "disturbance", change.date := decimal_date(max.disturbance.start.date)]
snv.changes[Classified.LT.Change == "step", change.date := decimal_date(breakpoint.date)]
#
# ---- make barplot ----

change.patterns <- snv.stats[ , .(num.genomes = .N, sum.abund = sum(sum.abund, na.rm = T), mean.abund = sum(mean.abund, na.rm = T), max.abund = sum(max.abund, na.rm = T)), by = .(phylum, Classified.LT.Change)]

tot.abunds <- change.patterns[ ,.(tot = sum(mean.abund)), by = .(phylum)]
tot.abunds <- tot.abunds[order(tot, decreasing = TRUE)]
tot.abunds[ ,order := 1:nrow(tot.abunds)]
change.patterns <- merge(change.patterns, tot.abunds, "phylum")
change.patterns <- change.patterns[ ,Classified.LT.Change := factor(Classified.LT.Change, levels = c("none", "disturbance","step","gradual"), ordered = TRUE)]
change.patterns <- change.patterns[order(order, Classified.LT.Change)]

my.mat <- dcast(data = change.patterns, formula = Classified.LT.Change ~ phylum, value.var = "mean.abund", fill = 0)
my.mat <- as.matrix(my.mat, rownames = TRUE)
col.index <- data.table("colname" = colnames(my.mat),"colorder" = 1:ncol(my.mat))
col.index <- merge(x = col.index, y = tot.abunds, by.x = "colname", by.y = "phylum")
col.index <- col.index[order(order),colorder]
my.mat <- my.mat[ ,col.index]

y.max <- max(colSums(my.mat))

par(mar = c(4,2.5,1,0))
bar.spots <- barplot(height = my.mat, legend.text = c("None","Disturbance/Resilience","Step Change","Gradual Change"), col = c("grey","gold2","red2","skyblue2"), axes = FALSE, ann = FALSE, names.arg = rep("",ncol(my.mat)))
text(x = bar.spots, y = -(y.max / 75), labels = sub(pattern = "p__", replacement = "", x = colnames(my.mat)), cex = 1, adj = 1, srt = 30, xpd = NA)
axis(side = 2, labels = F, line = -1)
axis(side = 2, lwd = 0, line = -1.25, las = 2)
mtext("Mean Abundance (%)", side = 2, line = 1.25)

# ---- make barplot with other ----

change.patterns <- snv.stats[ , .(num.genomes = .N, sum.abund = sum(sum.abund, na.rm = T), mean.abund = sum(mean.abund, na.rm = T), max.abund = sum(max.abund, na.rm = T)), by = .(phylum, Classified.LT.Change)]

tot.abunds <- change.patterns[ ,.(tot = sum(mean.abund)), by = .(phylum)]
tot.abunds <- tot.abunds[order(tot, decreasing = TRUE)]

tot.abunds[ , make.other := c(phylum[1:5], rep("Other",nrow(tot.abunds) - 5))]
change.patterns <- merge(change.patterns,tot.abunds,"phylum")
change.patterns <- change.patterns[ , .(num.genomes = .N, sum.abund = sum(sum.abund, na.rm = T), mean.abund = sum(mean.abund, na.rm = T), max.abund = sum(max.abund, na.rm = T)), by = .(make.other, Classified.LT.Change)]
tot.abunds <- tot.abunds[ ,.(tot = sum(tot)), by = .(make.other)]

tot.abunds[ ,order := 1:nrow(tot.abunds)]
change.patterns <- merge(change.patterns, tot.abunds, "make.other")
change.patterns <- change.patterns[ ,Classified.LT.Change := factor(Classified.LT.Change, levels = c("none", "disturbance","step","gradual"), ordered = TRUE)]
change.patterns <- change.patterns[order(order, Classified.LT.Change)]

my.mat <- dcast(data = change.patterns, formula = Classified.LT.Change ~ make.other, value.var = "mean.abund", fill = 0)
my.mat <- as.matrix(my.mat, rownames = TRUE)
col.index <- data.table("colname" = colnames(my.mat),"colorder" = 1:ncol(my.mat))
col.index <- merge(x = col.index, y = tot.abunds, by.x = "colname", by.y = "make.other")
col.index <- col.index[order(order),colorder]
my.mat <- my.mat[ ,col.index]

y.max <- max(colSums(my.mat))

par(mar = c(4,2.5,1,0))
bar.spots <- barplot(height = my.mat, legend.text = c("None","Disturbance/Resilience","Step Change","Gradual Change"), col = c("grey","gold2","red2","skyblue2"), axes = FALSE, ann = FALSE, names.arg = rep("",ncol(my.mat)))
text(x = bar.spots, y = -(y.max / 75), labels = sub(pattern = "p__", replacement = "", x = colnames(my.mat)), cex = 1, adj = 1, srt = 30, xpd = NA)
axis(side = 2, labels = F, line = -1)
axis(side = 2, lwd = 0, line = -1.25, las = 2)
mtext("Mean Abundance (%)", side = 2, line = 1.25)

# ---- make stripchart ----

# snv.changes <- snv.changes[ ,.(phylum, class, order, family, genus, species, mean.abund, med.abund, Classified.LT.Change, change.date)]
# 
# stripchart(x = change.date ~ genus, data = snv.changes[phylum == "p__Actinobacteriota"], vertical = TRUE, method = "jitter", pch = 21, jitter = .25, cex = )

snv.changes[ ,phylum := sub("p__","",phylum)]
snv.changes[ ,family := sub("f__","",family)]
snv.changes[ ,genus := sub("g__","",genus)]

acI <- snv.changes[family == "Nanopelagicaceae"]
acI[order(change.date)]

library(ggplot2)

ggplot(data = snv.changes, aes(y = change.date, x = phylum))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(), 
        axis.text.x.bottom = element_text(angle = 30))+
  geom_jitter(aes(color = Classified.LT.Change, size = mean.abund), alpha = .5)+
  scale_color_manual(values = c("gold2","red2"), labels = c("Disturbance/\nResilience", "Step Change"))+
  scale_y_continuous(breaks = 2000:2019)+
  guides(color = guide_legend(title = "Abrupt Change", override.aes = c(size = 5)))+
  guides(size = guide_legend(title = "Mean\nAbundance (%)"))

ggplot(data = snv.changes, aes(y = change.date, x = genus))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(), 
        axis.text.x.bottom = element_text(angle = 30))+
  geom_jitter(aes(color = Classified.LT.Change, size = mean.abund), alpha = .5)+
  scale_color_manual(values = c("gold2","red2"), labels = c("Disturbance/\nResilience", "Step Change"))+
  scale_y_continuous(breaks = 2000:2019)+
  guides(color = guide_legend(title = "Abrupt Change", override.aes = c(size = 5)))+
  guides(size = guide_legend(title = "Mean\nAbundance (%)"))

ggplot(data = snv.changes[phylum == "Actinobacteriota"], aes(y = change.date, x = genus))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(), 
        axis.text.x.bottom = element_text(angle = 30))+
  geom_jitter(aes(color = Classified.LT.Change, size = mean.abund), alpha = .5)+
  scale_color_manual(values = c("gold2","red2"), labels = c("Disturbance/\nResilience", "Step Change"))+
  scale_y_continuous(breaks = 2000:2019)+
  guides(color = guide_legend(title = "Abrupt Change", override.aes = c(size = 5)))+
  guides(size = guide_legend(title = "Mean\nAbundance (%)"))

ggplot(data = snv.changes, aes(y = change.date, x = family))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(), 
        axis.text.x.bottom = element_text(angle = 30))+
  geom_jitter(aes(color = Classified.LT.Change, size = mean.abund), alpha = .5)+
  scale_color_manual(values = c("gold2","red2"), labels = c("Disturbance/\nResilience", "Step Change"))+
  scale_y_continuous(breaks = 2000:2019)+
  guides(color = guide_legend(title = "Abrupt Change", override.aes = c(size = 5)))+
  guides(size = guide_legend(title = "Mean\nAbundance (%)"))

# ---- make as a beeswarm instead? ----


snv.changes[family == "Nanopelagicaceae", .(genome,genus,species,Classified.LT.Change,change.date)]
