# RRR
# Get SNV pres-abs tables
# goal, track when new SNVs arrive/become fixed
# but in this script just make the tables

# note: the /tmp user limit is 100GB on helheim
# and in each process, the SNV file forms a tmp file while the process is running
# so to not fill up the tmp folder and stop running, tell R to use a different tmp folder that doesn't have the same system limits in place
# but remember, this will take off a "safety" and you could more easily fill up disk space by accident. but anyway, 
# add to the .Renviron file a line that says:
# TMP = '/home/rrohwer/myRtempfiles/'
# now all R sessions will save tmp files here instead.
# nevermind that didn't work. Tell PARALLEL to use the different tmp directory:
# cat get_new_snvs.jobfile | parallel -j 60 --tmpdir ../myRtempfiles/ > termout_get_stats.txt 2>&1

# ---- set-up ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(lubridate))

userprefs <- commandArgs(trailingOnly = TRUE)
per.genome.snv.file <- userprefs[1]
genome.info.file <- userprefs[2]
num.threads <- as.numeric(userprefs[3])
output.folder <- userprefs[4]

# # local path testing
# per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-21_3300043464_group3_bin69_SNVs.tsv.gz" # B
# genome.info.file <- "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz"
# num.threads <- 1
# output.folder <- "data/2023-12-14_new_SNV_testing/"

# ---- filter to abund samples same as for multivariate SNV analyses & import data ----

# NOTE: that calc dist matrix script had breadth cutoff .5 instead of .7, but it shouldn't matter b/c med cov > 10x is way more stringent

my.genome <- sub(pattern = "^.*/", replacement = "", x = per.genome.snv.file)
my.genome <- sub(pattern = "_SNVs.*$", replacement = "", x = my.genome)

genome.info <- fread(file = genome.info.file, nThread = num.threads, colClasses = c("date" = "character"))
genome.info <- genome.info[genome == my.genome]

# Remove samples with coverage too low to reliably call SNVs
genome.info <- genome.info[ ,above.breadth.cutoff := (breadth / breadth_expected) >= .7]
genome.info <- genome.info[above.breadth.cutoff == TRUE & coverage_median > 10, ]

# Skip genomes that don't occur on at least 30 dates in at least 10 years
if (length(unique(genome.info$date)) < 30 | length(unique(genome.info$year)) < 10 ){
  quit(save = "no", status = 0)
}

cat("Processing", my.genome,"\n")

snv.instrain <- fread(file = per.genome.snv.file, nThread = num.threads, colClasses = c("date" = "character"))

snv.instrain <- snv.instrain[sample %in% unique(genome.info$sample)]

# ---- get total SNVs per day stats ----

snv.count <- snv.instrain[ , .(tot.SNVs = .N), by = .(sample, date, mutation_type)]
snv.count <- snv.count[ , .(tot.SNVs = mean(tot.SNVs)), by = .(date, mutation_type)]
tot.snv <- snv.count[ , .(tot.SNVs = sum(tot.SNVs)), by = .(date)]
tot.N <- snv.count[mutation_type == "N" , .(tot.N = sum(tot.SNVs)), by = .(date)]
tot.S <- snv.count[mutation_type == "S" , .(tot.S = sum(tot.SNVs)), by = .(date)]
snv.count <- merge(x = tot.snv, y = tot.N, by = "date")
snv.count <- merge(x = snv.count, y = tot.S, by = "date")

# ---- parse into pres/abs table for SNVs ----

snv.instrain[ ,is.pres := TRUE]

snv <- dcast(data = snv.instrain, formula = scaffold + position + gene + mutation_type ~ date, fun.aggregate = any, value.var = "is.pres", fill = FALSE)

snv[ ,snv.num := paste0("snv",1:nrow(snv))]
snv.key <- snv[ ,.(scaffold, position, snv.num, gene, mutation_type)]
snv[ ,`:=`(scaffold = NULL, position = NULL, gene = NULL, mutation_type = NULL)]
# setcolorder(x = snv, neworder = c(ncol(snv), 1:(ncol(snv) - 1)))

# ---- parse into new/recurring matrix ----

# "new" SNV appearances are when it steps from 0 to 1
# subtracting offset matrices can ID this

snv.mat <- as.matrix(x = snv, rownames = "snv.num")

snv.mat <- snv.mat[ ,-1] - snv.mat[ ,-ncol(snv.mat)]

snv.mat <- snv.mat == 1

# "recurring" SNVs are when it's new but the cumulative sum of times it was new is greater than 1
snv.mat.2 <- apply(X = snv.mat, MARGIN = 1, FUN = cumsum)
snv.mat.2 <- t(snv.mat.2)
# snv.mat[1:10,1:15]
# snv.mat.2[1:10,1:15]

my.mat <- matrix(data = NA, nrow = nrow(snv.mat), ncol = ncol(snv.mat), dimnames = list(row.names(snv.mat), colnames(snv.mat)))
my.mat[snv.mat == TRUE & snv.mat.2 == 1] <- "N"
my.mat[snv.mat == TRUE & snv.mat.2 != 1] <- "R"

# my.mat[1:10,1:15]

# ---- get number new/ number recurring per day ----

my.nonsyn <- my.mat[snv.key$mutation_type == "N", ]
my.syn  <- my.mat[snv.key$mutation_type == "S", ]

snv.stats <- data.table("Date" = colnames(my.mat), 
                        "New" = colSums(my.mat == "N", na.rm = T), 
                        "Recurring" = colSums(my.mat == "R", na.rm = T))
snv.stats[ ,`:=`(Either = New + Recurring,
                 New.Nonsyn = colSums(my.nonsyn == "N", na.rm = T),
                 Recurring.Nonsyn  = colSums(my.nonsyn == "R", na.rm = T))]
snv.stats[ ,`:=`(Either.Nonsyn = New.Nonsyn + Recurring.Nonsyn,
                 New.Syn = colSums(my.syn == "N", na.rm = T),
                 Recurring.Syn  = colSums(my.syn == "R", na.rm = T))]
snv.stats[ ,`:=`(Either.Syn = New.Syn + Recurring.Syn)]

# ---- combine stats and export ----

snv.stats <- merge(x = snv.stats, y = snv.count, by.x = "Date", by.y = "date", all = TRUE)

dates.key <- unique(snv.instrain[ ,.(date, year, yday, season, invasion)])
snv.stats <- merge(x = snv.stats, y = dates.key, by.x = "Date", by.y = "date", all = T) 

snv.stats[ ,genome := my.genome]

setcolorder(x = snv.stats, neworder = c("genome","Date","year","yday","season","invasion","New","Recurring","Either","New.Nonsyn","Recurring.Nonsyn","Either.Nonsyn","New.Syn","Recurring.Syn","Either.Syn","tot.SNVs","tot.N","tot.S"))
snv.stats <- snv.stats[order(parse_date_time(Date,"ymd"))] # this doesn't seem to be necessary, but just in case

fwrite(x = snv.stats, file = file.path(output.folder,paste0(my.genome,"_new_and_total_SNVs.tsv.gz")), sep = "\t", nThread = num.threads)

