# RRR

# my last script spit out a bunch of one-line tables. Need to simply combine them all into 1 table with all the included genomes

# OK this is actually works for the distance stats files too... just run it with different input. 
# Each genomes's distance stats tables are not that big and have a genome column in them already. 
# Need to simply combine them all into 1 table with all the included genomes
# that will be less unwieldy and the files will not be too large

userinput <- commandArgs(trailingOnly = TRUE)
folder.path <- userinput[1]
output.file <- userinput[2]
threads <- as.numeric(userinput[3])

# # test locally
# folder.path <- "data/2023-11-01_multidimensional_SNV_analysis/time_decay_example_plots"
# output.file <- "data/2023-11-01_multidimensional_SNV_analysis/time_decay_example_plots/example_combined_stats_file.tsv.gz"
# threads <- 1

# ---- set up ----

library(data.table)

# ---- import and combine ----

all.files <- list.files(path = folder.path, pattern = "*.tsv.*", full.names = T)

combo.list <- list()
for (f in all.files){
  my.tab <- fread(file = f, header = T, nThread = threads)
  combo.list <- c(combo.list, list(my.tab))
}

combo.table <- rbindlist(combo.list)

# ---- export ----

fwrite(x = combo.table, file = output.file, sep = "\t", nThread = threads)


# ---- run as: ----

# $ pwd
# /home/rrohwer/multivariate_SNV_analysis-med_cov_10

# $ Rscript 3b2_combine_time_decay_stats_files.R time_decay_stats_all_SNVs all_genomes_time_decay_stats-all_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R time_decay_stats_nonsynonymous_SNVs all_genomes_time_decay_stats-nonsynonymous_SNVs.tsv.gz 115

# $ Rscript 3b2_combine_time_decay_stats_files.R stats_all_SNVs/pairwise-season all_genomes_pairwise_season_distance_stats-all_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R stats_all_SNVs/pairwise-year all_genomes_pairwise_year_distance_stats-all_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R stats_all_SNVs/per-sample all_genomes_per_sample_distance_stats-all_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R stats_all_SNVs/per-season all_genomes_per_season_distance_stats-all_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R stats_all_SNVs/per-year all_genomes_per_year_distance_stats-all_SNVs.tsv.gz 115

# $ Rscript 3b2_combine_time_decay_stats_files.R stats_nonsynon_SNVs/pairwise-season all_genomes_pairwise_season_distance_stats-nonsynonymous_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R stats_nonsynon_SNVs/pairwise-year all_genomes_pairwise_year_distance_stats-nonsynonymous_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R stats_nonsynon_SNVs/per-sample all_genomes_per_sample_distance_stats-nonsynonymous_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R stats_nonsynon_SNVs/per-season all_genomes_per_season_distance_stats-nonsynonymous_SNVs.tsv.gz 115
# $ Rscript 3b2_combine_time_decay_stats_files.R stats_nonsynon_SNVs/per-year all_genomes_per_year_distance_stats-nonsynonymous_SNVs.tsv.gz 115




