# RRR
# jobfiles to run on helheim using GNU parallel


genomes <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
genomes <- genomes[genomes$winner, "bin.full.name"]

# for calc dist matrix and NMDS ----

per.genome.snv.file <- paste0("per-genome_SNVs/",genomes,"_SNVs.tsv.gz")
combo.genome.info.file <- "../yggshare/current_projects/TYMEFLIES/tymeflies/processinstrain96/genome_info_combined.tsv.gz"
output.dist.folder.all <- "multivariate_SNV_analysis-med_cov_10/distance_matrices_all_SNVs"
output.dist.folder.N <- "multivariate_SNV_analysis-med_cov_10/distance_matrices_nonsynon_SNVs"
output.nmds.object.folder.all <- "multivariate_SNV_analysis-med_cov_10/NMDS_objects_all_SNVs"
output.nmds.object.folder.N <- "multivariate_SNV_analysis-med_cov_10/NMDS_objects_nonsynon_SNVs"
num.threads <- 1

step1.jobfile <- paste0("Rscript 1-calculate_distance_matrix_and_NMDS.R ",
                        per.genome.snv.file," ",
                        combo.genome.info.file," ",
                        output.dist.folder.all," ",
                        output.dist.folder.N," ",
                        output.nmds.object.folder.all," ",
                        output.nmds.object.folder.N," ",
                        1)

write.table(x = step1.jobfile, file = "server-scripts/generated/2023-11-01_multivariate_SNV_analysis/calc_distance_and_nmds.jobfile", quote = F, row.names = F, col.names = F)


# for calc multivariate stats - all SNVs ----

dist.file <- paste0("multivariate_SNV_analysis-med_cov_10/distance_matrices_all_SNVs/",genomes,"_all_SNV_euclidean_distance_matrix.rds")
sample.key <- "sample_key.tsv"
output.per.sample.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_all_SNVs/per-sample/",genomes,"_all_SNV_per_sample_stats.tsv.gz")
output.per.year.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_all_SNVs/per-year/",genomes,"_all_SNV_per_year_stats.tsv.gz")
output.pairwise.year.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_all_SNVs/pairwise-year/",genomes,"_all_SNV_pairwise_year_stats.tsv.gz")
output.per.season.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_all_SNVs/per-season/",genomes,"_all_SNV_per_season_stats.tsv.gz")
output.pairwise.season.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_all_SNVs/pairwise-season/",genomes,"_all_SNV_pairwise_season_stats.tsv.gz")
threads <- 1

step.2b.jobfile <- paste("Rscript 2b-get_distance_stats.R",
                         dist.file,
                         sample.key,
                         output.per.sample.stats,
                         output.per.year.stats,
                         output.pairwise.year.stats,
                         output.per.season.stats,
                         output.pairwise.season.stats,
                         1)

write.table(x = step.2b.jobfile, file = "server-scripts/generated/2023-11-01_multivariate_SNV_analysis/calc_distance_stats-all_SNVs.jobfile", quote = F, row.names = F, col.names = F)

# for calc multivariate stats - nonsynonymous SNVs ----

dist.file <- paste0("multivariate_SNV_analysis-med_cov_10/distance_matrices_nonsynon_SNVs/",genomes,"_nonsynonymous_SNV_euclidean_distance_matrix.rds")
sample.key <- "sample_key.tsv"
output.per.sample.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_nonsynon_SNVs/per-sample/",genomes,"_nonsynonymous_SNV_per_sample_stats.tsv.gz")
output.per.year.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_nonsynon_SNVs/per-year/",genomes,"_nonsynonymous_SNV_per_year_stats.tsv.gz")
output.pairwise.year.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_nonsynon_SNVs/pairwise-year/",genomes,"_nonsynonymous_SNV_pairwise_year_stats.tsv.gz")
output.per.season.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_nonsynon_SNVs/per-season/",genomes,"_nonsynonymous_SNV_per_season_stats.tsv.gz")
output.pairwise.season.stats <- paste0("multivariate_SNV_analysis-med_cov_10/stats_nonsynon_SNVs/pairwise-season/",genomes,"_nonsynonymous_SNV_pairwise_season_stats.tsv.gz")
threads <- 1

step.2b.jobfile.nonsyn <- paste("Rscript 2b-get_distance_stats.R",
                         dist.file,
                         sample.key,
                         output.per.sample.stats,
                         output.per.year.stats,
                         output.pairwise.year.stats,
                         output.per.season.stats,
                         output.pairwise.season.stats,
                         1)

write.table(x = step.2b.jobfile.nonsyn, file = "server-scripts/generated/2023-11-01_multivariate_SNV_analysis/calc_distance_stats-nonsynonymous_SNVs.jobfile", quote = F, row.names = F, col.names = F)

# to combine the all-genome stats files, just concatenate them. can I do it without unzipping? like
# cat stats_all_SNVs/per-sample/* > stats_genomes_combined/all_genomes_all_SNVs_per_sample_distance_stats.tsv.gz
# seems like it works just fine

# for plot time decays - all SNVs ----

dist.file <- paste0("distance_matrices_all_SNVs/",genomes,"_all_SNV_euclidean_distance_matrix.rds")
sample.key <- "sample_key.tsv"
tax.file <- "drep_results_all_genomes_0.96.rds"
output.plot.folder <- "time_decay_plots_all_SNVs"
output.stats.folder <- "time_decay_stats_all_SNVs"
threads <- 1

step.3b.jobfile.all <- paste("Rscript 3b-make_time_decay_plots.R",
                             dist.file,
                             sample.key,
                             tax.file,
                             output.plot.folder,
                             output.stats.folder,
                             threads)

write.table(x = step.3b.jobfile.all, file = "server-scripts/generated/2023-11-01_multivariate_SNV_analysis/do_all_SNV_time_decays.jobfile", quote = F, row.names = F, col.names = F)

# for plot time decays - nonsynonymous SNVs ----

dist.file <- paste0("distance_matrices_nonsynon_SNVs/",genomes,"_nonsynonymous_SNV_euclidean_distance_matrix.rds")
sample.key <- "sample_key.tsv"
tax.file <- "drep_results_all_genomes_0.96.rds"
output.plot.folder <- "time_decay_plots_nonsynonymous_SNVs"
output.stats.folder <- "time_decay_stats_nonsynonymous_SNVs"
threads <- 1

step.3b.jobfile.nonsyn <- paste("Rscript 3b-make_time_decay_plots.R",
                             dist.file,
                             sample.key,
                             tax.file,
                             output.plot.folder,
                             output.stats.folder,
                             threads)

write.table(x = step.3b.jobfile.nonsyn, file = "server-scripts/generated/2023-11-01_multivariate_SNV_analysis/do_nonsyn_SNV_time_decays.jobfile", quote = F, row.names = F, col.names = F)

# for make NMDS plots- all SNVs ----

nmds.file <- paste0("NMDS_objects_all_SNVs/",genomes,"_all_SNV_euclidean_distance_NMDS_object.rds")
sample.key.file <- "sample_key.tsv"
selection.summary.file <- paste0("../../yggshare/current_projects/TYMEFLIES/tymeflies/selected_genes_analysis/per-genome_selection_summaries/",genomes,"_selection_summary.tsv.gz")
genome.file <- "../../yggshare/current_projects/TYMEFLIES/tymeflies/processinstrain96/genome_info_combined.tsv.gz"
tax.file <- "drep_results_all_genomes_0.96.rds"
output.folder.tables <- "NMDS_tables_all_SNVs"
output.folder.plots <- "NMDS_plots_all_SNVs"
threads <-1

step.2a.jobfile.all <- paste("Rscript 2a-make_NMDS_plots.R",
                             nmds.file,
                             sample.key.file,
                             selection.summary.file,
                             genome.file,
                             tax.file,
                             output.folder.tables,
                             output.folder.plots,
                             threads)

write.table(x = step.2a.jobfile.all, file = "server-scripts/generated/2023-11-01_multivariate_SNV_analysis/make_NMDS_plots-all_SNVs.jobfile", quote = F, row.names = F, col.names = F)


# for make NMDS plots- nonsynonymous SNVs ----

nmds.file <- paste0("NMDS_objects_nonsynon_SNVs/",genomes,"_nonsynonymous_SNV_euclidean_distance_NMDS_object.rds")
sample.key.file <- "sample_key.tsv"
selection.summary.file <- paste0("../../yggshare/current_projects/TYMEFLIES/tymeflies/selected_genes_analysis/per-genome_selection_summaries/",genomes,"_selection_summary.tsv.gz")
genome.file <- "../../yggshare/current_projects/TYMEFLIES/tymeflies/processinstrain96/genome_info_combined.tsv.gz"
tax.file <- "drep_results_all_genomes_0.96.rds"
output.folder.tables <- "NMDS_tables_nonsynonymous_SNVs"
output.folder.plots <- "NMDS_plots_nonsynonymous_SNVs"
threads <-1

step.2a.jobfile.non <- paste("Rscript 2a-make_NMDS_plots.R",
                             nmds.file,
                             sample.key.file,
                             selection.summary.file,
                             genome.file,
                             tax.file,
                             output.folder.tables,
                             output.folder.plots,
                             threads)

write.table(x = step.2a.jobfile.non, file = "server-scripts/generated/2023-11-01_multivariate_SNV_analysis/make_NMDS_plots-nonsyn_SNVs.jobfile", quote = F, row.names = F, col.names = F)





