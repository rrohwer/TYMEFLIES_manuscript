# RRR

# need a folder with just the WINNER bin fasta in it
# this is for singleM input
# probably good to have handy anyway
# before I always used them concatenated


bin.stats <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
colnames(bin.stats)

# ---- make script ----

bin.stats <- bin.stats[bin.stats$winner, ]

# will run from outside the metabat2_bins_renamed folder
bin.locs <- paste("metabat2_bins_renamed",bin.stats$tymeflies.name, bin.stats$bin.filename, sep = "/")
new.locs <- paste0("representative_genome_fastas/",bin.stats$bin.filename)


bash.copy.script <- c("#!/bin/bash",
                        paste0("cp ",bin.locs," ",new.locs))

# ---- export script ----

write.table(x = bash.copy.script, file = "server-scripts/generated/2023-12-18_singleM/copy_bin_fastas.sh",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
