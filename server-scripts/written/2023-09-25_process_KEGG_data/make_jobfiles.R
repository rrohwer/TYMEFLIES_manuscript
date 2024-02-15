# RRR
# parallel jobfiles for keg R script processing

bin.names <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
bin.names <- bin.names$bin.full.name[bin.names$winner]


# for adding pathway and module info:
kofam.files <- paste0("../yggshare/current_projects/TYMEFLIES/tymeflies/kegg_annotations/kofamscan_kegg_annotations/",bin.names,".kofamscan.tsv.gz")
output.one.line.per.gene <- paste0("kegg_annotations_with_pathways-one_line_per_gene/",bin.names,".ko-pathways_genes-unique.tsv.gz")
output.one.line.per.ann <- paste0("kegg_annotations_with_pathways/",bin.names,".ko-pathways_anns-separate.tsv.gz")

add.path.jobfile <- paste0("Rscript add_pathway_info_to_KO_annotations.R ",
                           "KEGG_lists/ko_pathway_links.txt ",
                           "KEGG_lists/ko_module_links.txt ",
                           "KEGG_lists/module_pathway_links.txt ",
                           "KEGG_lists/pathway_list.txt ",
                           "KEGG_lists/ko_list.txt ",
                           "KEGG_lists/module_list.txt ",
                           kofam.files," ",
                           output.one.line.per.gene," ",
                           output.one.line.per.ann," ",
                           "1")

write.table(x = add.path.jobfile, file = "server-scripts/generated/2023-09-25_process_KEGG_data/add_pathway_to_ko.jobfile", quote = F, row.names = F, col.names = F)


# for finding consistently selected genes:
gene.info.files <- paste0("../yggshare/current_projects/TYMEFLIES/tymeflies/processinstrain96/per-genome_gene_info_MK/",bin.names,"_gene_info_combined_MK.tsv.gz")
annotation.files <- paste0("../yggshare/current_projects/TYMEFLIES/tymeflies/kegg_annotations/kegg_annotations_with_pathways-one_line_per_gene/",bin.names,".ko-pathways_genes-unique.tsv.gz")
coverm.file <- "coverM_drep96_minANI93.txt.gz"
genome.file <- "../yggshare/current_projects/TYMEFLIES/tymeflies/processinstrain96/genome_info_combined.tsv.gz"
output.key <- paste0("consistent_genes_analysis/key_genes/",bin.names,".consistently_selected_genes.tsv.gz")
output.daily <- paste0("consistent_genes_analysis/key_genes-by_sample/",bin.names,".consistently_selected_genes-by_sample.tsv.gz")

key.genes.jobfile <- paste("Rscript find_consistently_selected_genes.R",
                            gene.info.files,
                            annotation.files,
                            coverm.file,
                            genome.file,
                            bin.names,
                            output.key,
                            output.daily,
                            1, sep = " ")

write.table(x = key.genes.jobfile, file = "server-scripts/generated/2023-09-25_process_KEGG_data/find_key_genes.jobfile", quote = F, row.names = F, col.names = F)


# for making plots:
genes.files <- paste0("consistent_genes_analysis/key_genes-by_sample/",bin.names,".consistently_selected_genes-by_sample.tsv.gz") #output.daily from above is still in home dir

plotting.jobfile <- paste0("Rscript baseR_heatmaps_of_selected_genes.R ",
                           genes.files, 
                           " 1 ",
                           "consistently_selected_genes_heatmaps")

write.table(x = plotting.jobfile, file = "server-scripts/generated/2023-09-25_process_KEGG_data/plot_key_genes.jobfile", quote = F, row.names = F, col.names = F)
