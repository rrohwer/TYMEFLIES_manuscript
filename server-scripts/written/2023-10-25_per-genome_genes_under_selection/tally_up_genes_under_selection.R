# RRR
# scale up my corrected genes under selection plots to make them for all genomes

# ---- set-up ----

library(data.table)
library(lubridate)

userinput <- commandArgs(trailingOnly = TRUE)
genes <- userinput[1]
output <- userinput[2]
num.threads <- as.numeric(userinput[3])

# # local paths for testing
# genes <- "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz"
# output <- "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_selection_summary.tsv.gz"
# num.threads <- 1

# ---- parse data ----

genome <- sub(pattern = "^.*/", replacement = "", x = genes)
genome <- sub(pattern = "_gene_info_combined_MK.tsv.gz", replacement = "", x = genome)

genes <- fread(file = genes, nThread = num.threads)

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,sample := sub(pattern = "\\.IS_gene_info\\.tsv", replacement = "", x = sample)]
genes[ ,sample := sub(pattern = ".gz", replacement = "", x = sample)]

pos <- genes[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman), .(.N), by = .(sample, date, season, year)]
neg <- genes[MK.p.val <= .05 & mcdonald.kreitman > 1 & is.finite(mcdonald.kreitman), .(.N), by = .(sample, date, season, year)] 

sel.summary <- merge(x = pos, y = neg, by = c("sample", "date", "season", "year"), all = T, suffixes = c("pos","neg"))
sel.summary[is.na(Npos), Npos := 0]
sel.summary[is.na(Nneg), Nneg := 0]

sel.summary$genome <- genome

setcolorder(x = sel.summary, neworder = c("genome","sample","date","season","year","Npos","Nneg"))

sel.summary[ ,date := as.character(date)]

# ---- save by-genome files ----

fwrite(x = sel.summary, file = output, nThread = num.threads)
