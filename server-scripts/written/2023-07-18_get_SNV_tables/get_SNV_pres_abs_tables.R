# RRR
# Get SNV pres-abs tables
# goal, track when new SNVs arrive/become fixed
# but in this script just make the tables

# ---- set-up ----

library(data.table)

# local path testing
# cat("\n\nComment out your local paths!\n\n")
# per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-21_3300043464_group3_bin69_SNVs.tsv.gz" # B
# per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2016-07-20_3300033996_group7_bin32_SNVs.tsv.gz" # C
# per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-04_3300044729_group3_bin142_SNVs.tsv.gz" # A
# num.threads <- 1
# output.pres.abs.folder <- "data/2023-07-16_SNVs_example_files/generated_pres_abs"

# on server
userprefs <- commandArgs(trailingOnly = TRUE)
per.genome.snv.file <- userprefs[1]
num.threads <- as.numeric(userprefs[2])
output.pres.abs.folder <- userprefs[3] # make before running script!
output.snv.stats.folder <- userprefs[4] # make before running script!

cat("working on", per.genome.snv.file, "\n")

# ---- parse into pres/abs table for SNVs ----

x <- "character"
names(x) <- "date"
snv.instrain <- fread(file = per.genome.snv.file, nThread = num.threads, colClasses = x)
my.genome <- snv.instrain[1, genome]
snv.instrain[ ,is.pres := TRUE]

snv <- dcast(data = snv.instrain, formula = scaffold + position + gene + mutation_type ~ date, fun.aggregate = any, value.var = "is.pres", fill = FALSE)

snv[ ,snv.num := paste0("snv",1:nrow(snv))]
snv.key <- snv[ ,.(scaffold, position, snv.num, gene, mutation_type)]
snv[ ,`:=`(scaffold = NULL, position = NULL, gene = NULL, mutation_type = NULL)]
setcolorder(x = snv, neworder = c(ncol(snv), 1:(ncol(snv) - 1)))

# ---- export pres/abs table ----

fwrite(x = snv, file = file.path(output.pres.abs.folder, paste0(my.genome, "_SNV_pres_abs.tsv.gz")), compress = "gzip", sep = "\t", nThread = num.threads)
fwrite(x = snv.key, file = file.path(output.pres.abs.folder, paste0(my.genome, "_SNV_key.tsv.gz")), compress = "gzip", sep = "\t", nThread = num.threads)

# ---- parse into new/recurring matrix ----

# "new" SNV appearances are when it steps from 0 to 1
# subtracting offset matrices can ID this

snv.mat <- as.matrix(x = snv, rownames = "snv.num")

snv.mat <- snv.mat[ ,-1] - snv.mat[ ,-ncol(snv.mat)]

snv.mat <- snv.mat == 1

# "recurring" SNVs are when it's new but the cumulative sum of times it was new is greater than 1
snv.mat.2 <- apply(X = snv.mat, MARGIN = 1, FUN = cumsum)
snv.mat.2 <- t(snv.mat.2)

my.mat <- matrix(data = NA, nrow = nrow(snv.mat), ncol = ncol(snv.mat), dimnames = list(row.names(snv.mat), colnames(snv.mat)))
my.mat[snv.mat == TRUE & snv.mat.2 == 1] <- "N"
my.mat[snv.mat == TRUE & snv.mat.2 != 1] <- "R"

# ---- get number new/ number recurring per day ----

my.nonsyn <- my.mat[snv.key$mutation_type == "N", ]
my.syn  <- my.mat[snv.key$mutation_type == "S", ]

snv.stats <- data.table("Date" = colnames(my.mat), 
                        "Total.SNVs" = colSums(snv[ , -c(1,2)]), # snv stats missing day 1
                        "New" = colSums(my.mat == "N", na.rm = T), 
                        "Recurring" = colSums(my.mat == "R", na.rm = T))
snv.stats[ ,`:=`(Either = New + Recurring,
                 New.Nonsyn = colSums(my.nonsyn == "N", na.rm = T),
                 Recurring.Nonsyn  = colSums(my.nonsyn == "R", na.rm = T))]
snv.stats[ ,`:=`(Either.Nonsyn = New.Nonsyn + Recurring.Nonsyn,
                 New.Syn = colSums(my.syn == "N", na.rm = T),
                 Recurring.Syn  = colSums(my.syn == "R", na.rm = T))]
snv.stats[ ,`:=`(Either.Syn = New.Syn + Recurring.Syn)]


dates.key <- unique(snv.instrain[ ,.(date, year, yday, season, invasion)])

snv.stats <- merge(x = snv.stats, y = dates.key, by.x = "Date", by.y = "date", all.x = T, all.y = F)

# ---- export snv stats table ----

fwrite(x = snv.stats, file = file.path(output.snv.stats.folder, paste0(my.genome, "SNV_stats.tsv.gz")), compress = "gzip", sep = "\t", nThread = num.threads)

