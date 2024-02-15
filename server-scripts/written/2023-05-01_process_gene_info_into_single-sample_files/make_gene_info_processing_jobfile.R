# RRR

# example syntax:
# Rscript process_one-genome_gene_info.R ME2013-02-02s1D0_3300042325_group4_bin39 ../../yggshare/current_projects/TYMEFLIES/tymeflies/runinstrain96 drep_results_all_genomes_0.96.rds limony_seasons.rds 65

all.bins <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
all.bins <- all.bins[all.bins$winner, "bin.full.name", drop = T]

jobfile <- paste0("Rscript process_one-genome_gene_info.R ",all.bins," ../../yggshare/current_projects/TYMEFLIES/tymeflies/runinstrain96 drep_results_all_genomes_0.96.rds limony_seasons.rds 1")

write.table(x = jobfile, file = "server-scripts/generated/2023-05-01_process_gene_info_into_single-sample_files/make_gene_info_files.jobfile",
            quote = F, row.names = F, col.names = F)
