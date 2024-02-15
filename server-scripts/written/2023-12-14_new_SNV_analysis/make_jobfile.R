# RRR
# jobfiles to run on helheim using GNU parallel


genomes <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
genomes <- genomes[genomes$winner, "bin.full.name"]

# ---- get new-SNV stats ----

per.genome.snv.file <- paste0("../per-genome_SNVs/",genomes,"_SNVs.tsv.gz")
genome.info.file <- "../../yggshare/current_projects/TYMEFLIES/tymeflies/processinstrain96/genome_info_combined.tsv.gz"
num.threads <- 1
output.folder <- "new_SNV_stats"

get_new_snvs.jobfile <- paste("Rscript tally_new_SNVs_per_sample.R",
                              per.genome.snv.file,
                              genome.info.file,
                              num.threads,
                              output.folder)

write.table(x = get_new_snvs.jobfile, file = "server-scripts/generated/2023-12-14_new_SNV_analysis/get_new_snvs.jobfile", quote = F, row.names = F, col.names = F)

# ---- make plots ----

per.genome.stats.file <- paste0("new_SNV_stats/",genomes,"_new_and_total_SNVs.tsv.gz")
tax.file <- "drep_results_all_genomes_0.96.rds"
num.threads <- 1
output.folder <- "new_SNVs_plots"

plot_new_snvs.jobfile <- paste("Rscript plot_new_SNVs.R",
                               per.genome.stats.file,
                               tax.file,
                               num.threads,
                               output.folder)

write.table(x = plot_new_snvs.jobfile, file = "server-scripts/generated/2023-12-14_new_SNV_analysis/plot_new_snvs.jobfile", quote = F, row.names = F, col.names = F)
