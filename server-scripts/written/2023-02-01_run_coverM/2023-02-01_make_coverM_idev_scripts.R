# RRR
# make the coverM command that includes all the bam files
# set it up to run on the server, not sure how long it will take with all the files
# but was super fast with just 3 bams

# this was my attempt at making it for idev. Ended up taking too long so did not use this.

# ---- set up ----

sample.names <- readRDS(file = "data/2022-11-15_bin_stats/all_bins.rds")
sample.names <- unique(sample.names$tymeflies.name)

# ---- make scripts ----

all.bams <- paste0(sample.names,".bam")
all.bams <- paste0("../maptodrep/", all.bams)
all.bams <- paste0(all.bams, collapse = " ")

coverM.idev <- paste0("coverm genome --bam-files ",all.bams,
                      " --genome-definition coverM_genome_reference.txt ",
                      "--methods mean ",
                      "--contig-end-exclusion 5000 ",
                      "--output-file coverM_output_mean_5000.txt ",
                      "--output-format dense ",
                      "--threads 128 ",
                      "&> coverm_total_mapped_output_5000.txt",
                      " &") # save terminal output


# ---- export ----

write.table(x = coverM.idev, file = "server-scripts/generated/2023-02-01_run_coverM/idev.mean.5000.txt", 
            quote = F, row.names = F, col.names = F)
