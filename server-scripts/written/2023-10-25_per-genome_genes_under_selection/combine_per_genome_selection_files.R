# RRR

# this script is very simple
# I left the selection summaries as per-genome files so that I can easily make plots from them in parallel
# But I need a combined file of all the genome info to make my summary plots
# this combined file should not be too large, but will do on the server b/c that's where the per-genome files are

# ---- set-up ----

library(data.table)

userinput <- commandArgs(trailingOnly = TRUE)
input.folder <- userinput[1]
output.file <- userinput[2]
num.threads <- as.numeric(userinput[3])

# # local test-paths
# input.folder <- "data/2023-10-26_selection_per_sample_summaries/example_data"
# output.file <- "data/2023-10-26_selection_per_sample_summaries/example_combo_file.tsv.gz"
# num.threads <- 1


# ---- combine files into one ----

my.files <- list.files(path = input.folder, pattern = "^ME")

combo.list <- list()
for (f in my.files){
  temp <- fread(input = file.path(input.folder,f), nThread = num.threads, colClasses = c("character","character","character","character","numeric","numeric","numeric"))
  combo.list <- c(combo.list, list(temp))
}

combo.table <- rbindlist(l = combo.list)

# ---- save output file ----

fwrite(x = combo.table, file = output.file, nThread = num.threads)

