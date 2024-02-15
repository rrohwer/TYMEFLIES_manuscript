# RRR
# find the drep winners amongst all the bins
# concatenate them into a single fna file
# export a shell script to do this


bin.stats <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
colnames(bin.stats)

# ---- make script ----

bin.stats <- bin.stats[bin.stats$winner, ]

# will run from inside the metabat2_bins_renamed folder
bin.locs <- paste(bin.stats$tymeflies.name, bin.stats$bin.filename, sep = "/")

bin.locs

bash.concat.script <- c("#!/bin/bash",
                        paste0("cat ",bin.locs[1]," > drep_96_winners_concat.fna"),
                        paste0("cat ",bin.locs[-1]," >> drep_96_winners_concat.fna"))
bash.concat.script[1:10]

# ---- export script ----

write.table(x = bash.concat.script, file = "server-scripts/generated/2023-03-02-22_map_to_dRep_96_winners//1-concatenate_drep_winners.sh",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
