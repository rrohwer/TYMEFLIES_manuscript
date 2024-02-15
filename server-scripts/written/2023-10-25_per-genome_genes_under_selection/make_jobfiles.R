# RRR

genomes <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
genomes <- genomes[genomes$winner, "bin.full.name"]


# tally up # genes under selection per sample

genes.files <- paste0("per-genome_gene_info_MK/",genomes,"_gene_info_combined_MK.tsv.gz")
output.files <- paste0("per-genome_selection_summaries/",genomes,"_selection_summary.tsv.gz")

tally.jobfile <- paste0("Rscript tally_up_genes_under_selection.R ",genes.files," ",output.files," 1")

write.table(x = tally.jobfile, file = "server-scripts/generated/2023-10-25_per-genome_genes_under_selection/tally_up_genes_under_selection.jobfile", quote = F, row.names = F, col.names = F)
