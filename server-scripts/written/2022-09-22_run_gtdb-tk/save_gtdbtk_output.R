# RRR
# I ran checkM2 per folder, so need to move all the folders off of scratch
# Also need to zip the protein_files folders before trying to transfer them

bins <- readRDS("data/2022-09-14_list_of_all_bins/bin_key.rds")

sample.names <- unique(bins$sample)


# ---- save contents of the gtdb-tk error logs to check if any had errors (they didn't) ----

check.for.errors <- c("#!/bin/bash",
                      "echo gtdbtk error outputs > all_errors.txt",
                      "")

for (s in sample.names){
  check.for.errors <- c(check.for.errors,
                        "",
                        paste0("cd ", s),
                        paste0("echo ",s," >> ../all_errors.txt"),
                        paste0("cat gtdbtk.warnings.log >> ../all_errors.txt"),
                        "cd ../")
  
}

write.table(x = check.for.errors, file = "server-scripts/generated/2022-09-22_run_gtdb-tk/1-check_for_gtdbtk_errors.sh", append = F, quote = F, sep = "\n", row.names = F, col.names = F)

# ----

gather.output.files <- c("#!/bin/bash")

for (s in sample.names){
  gather.output.files <- c(gather.output.files,
                           "",
                           paste0("cd ",s),
                           paste0("cp gtdbtk.bac120.summary.tsv ",s,".tsv"),
                           paste0("mv ",s,".tsv ../gtdbtk_bac120_summary_files/"),
                           "cd ../")
}

write.table(x = gather.output.files, file = "server-scripts/generated/2022-09-22_run_gtdb-tk/2-put_all_output_summaries_into_single_folder.sh", append = F, quote = F, sep = "\n", row.names = F, col.names = F)