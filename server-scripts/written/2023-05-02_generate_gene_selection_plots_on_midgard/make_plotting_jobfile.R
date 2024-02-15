# RRR

# make the jobfiles for generating data tables and then also making the plots

all.bins <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
all.bins <- all.bins[all.bins$winner, "bin.full.name", drop = T]

# data syntax:
# Rscript find_selected_Cog_categories.R genes cogdefs anns genomes threads data.folder 

data.jobfile <- paste0("Rscript find_selected_COG_categories.R ",
                       "../../yggshare/current_projects/TYMEFLIES/tymeflies/processinstrain96/per-genome_gene_info/",all.bins,"_gene_info_combined.tsv.gz ",
                       "COG_descriptions.csv ",
                       "../../yggshare/current_projects/TYMEFLIES/tymeflies/eggnog_annotations/",all.bins,".emapper.annotations.gz ",
                       "drep_results_all_genomes_0.96.rds ",
                       "1 ",
                       "./")

write.table(x = data.jobfile, file = "server-scripts/generated/2023-05-02_generate_gene_selection_plots_on_midgard/add_MK_selection_stats_to_per-genome_gene_info_files.jobfile",
            quote = F, row.names = F, col.names = F)


# plot syntax:

#  Rscript plot_selected_COG_categories.R genes anns threads plot.folder

jobfile <- paste0("Rscript plot_selected_COG_categories.R ",
                  "../per-genome_gene_info_MK/",all.bins,"_gene_info_combined_MK.tsv.gz ",
                  "../../yggshare/current_projects/TYMEFLIES/tymeflies/eggnog_annotations/",all.bins,".emapper.annotations.gz ",
                  "1 ",
                  "./")

write.table(x = jobfile, file = "server-scripts/generated/2023-05-02_generate_gene_selection_plots_on_midgard/make_all_per-genome_selected_COGs_plots.jobfile",
            quote = F, row.names = F, col.names = F)


# ran with 
# cat make_all_per-genome_selected_COGs_plots.jobfile | parallel -j 30
# started: 4:55pm
