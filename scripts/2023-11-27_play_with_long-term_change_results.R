# RRR

# Things ran on the server. time decays and such. May not be perfect, may need to troubleshoot the breakpoints more.
# but first let's take a look at the output. 

library(data.table)
library(lubridate)
library(ggplot2)

snv.stats <- fread(file = "data/2023-11-01_multidimensional_SNV_analysis/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")
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

summary(tax$times.present)
summary(tax$times.present.med.cov.10)
summary(tax$mean.abund)

snv.stats <- merge(x = snv.stats, y = tax, by.x = "genome", by.y = "Genome")

snv.stats[ ,`:=`(max.disturbance.start.date = parse_date_time(max.disturbance.start.date, orders = "ymd"),
                breakpoint.date = parse_date_time(breakpoint.date, orders = "ymd"))]

# ---- classify into long-term change patterns ----

hist(snv.stats$interp.slope, breaks = 100)
abline(v=-2, col = "red")
# so if see a lot of -2 slopes, then should discount slopes less than positive 2 as well.
snv.stats$gradual.change <- FALSE
snv.stats[interp.slope > 2 & interp.adj.R2 > .5, gradual.change := TRUE]

snv.stats$step.change <- FALSE
snv.stats[!is.na(step.change.loc) &
            !is.na(breakpoint.loc) &
            step.change.length > 3.5 &
            gradual.change == FALSE, step.change := TRUE]
# this is calling too many right now. need to change things in the 3b script

snv.stats$disturbance.resilience <- FALSE
snv.stats[gradual.change == FALSE &
            step.change == FALSE &
            total.disturbances > 0, disturbance.resilience := TRUE]

snv.stats$lt.pattern <- "Stable"
snv.stats[gradual.change == TRUE, lt.pattern := "Gradual change"]
snv.stats[step.change == TRUE, lt.pattern := "Step change"]
snv.stats[disturbance.resilience == TRUE, lt.pattern := "Disturbance/Resilience"]
snv.stats[ ,lt.pattern <- factor(lt.pattern, levels = c("Gradual change", "Step change", "Disturbance/Resilience","Stable"), ordered = T)]

# pull out the year of disturbance or step changes

snv.changes <- snv.stats[disturbance.resilience == TRUE | step.change == TRUE]
snv.changes[disturbance.resilience == TRUE, change.year := year(max.disturbance.start.date)]
snv.changes[step.change == TRUE, change.year := year(breakpoint.date)]

# ---- make some quick plots ----

my.data <- snv.stats
my.data <- snv.stats[times.present.med.cov.10 > 20, ]


ggplot(data = my.data, aes(x = phylum, y = mean.abund, fill = lt.pattern))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_bar(stat = "identity")

# this is sooo ugly

ggplot(data = snv.changes, aes(y = change.year, x = phylum))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_jitter(aes(color = lt.pattern, size = mean.abund), alpha = .5)+
  geom_abline(intercept = 2012, slope = 0, color = adjustcolor("grey",.5))+
  geom_abline(intercept = 2014, slope = 0, color = adjustcolor("grey",.5))

ggplot(data = snv.changes[phylum == "p__Actinobacteriota"], aes(y = change.year, x = family))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_jitter(aes(color = lt.pattern, size = mean.abund), alpha = .5)+
  geom_abline(intercept = 2012, slope = 0, color = adjustcolor("grey",.5))+
  geom_abline(intercept = 2014, slope = 0, color = adjustcolor("grey",.5))











