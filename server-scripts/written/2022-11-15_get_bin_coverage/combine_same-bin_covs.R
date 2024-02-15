# RRR

# mapping 471 samples to each bin results in 471 cov files per bin.
# This script combines them into a single cov file, using a tall format

library(data.table)

userprefs <- commandArgs(trailingOnly = TRUE)
bin.folder <- userprefs[1]
results.folder <- userprefs[2]

cat("UNCOMMENT YOUR TROUBLESHOOTING LOCAL PATHS")
bin.folder <- "data/2022-06-02_cov_files_for_ME2000-03-15pf_3300044539/covfiles/"
results.folder

my.covs <- list.files(path = bin.folder, pattern = "\\.cov", full.names = F)

my.samples <- sub(pattern = "\\.cov", replacement = "", x = my.covs)

combo.list <- NULL
for (c in 1:length(my.covs)){
  combo.list[[c]] <- fread(file = file.path("data/2022-06-02_cov_files_for_ME2000-03-15pf_3300044539/covfiles/",my.covs[c]), header = T, sep = "\t")
  combo.list[[c]]$Sample <- my.samples[c]
}

combo.table <- do.call(what = rbind, args = combo.list)


str(combo.table)
colnames(combo.table)[1] <- "Contig"

# does data.table have a fast write too?