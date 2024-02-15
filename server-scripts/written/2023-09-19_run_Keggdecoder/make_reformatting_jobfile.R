# RRR

key <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
key <- key[key$winner,"bin.full.name"]
key.in <- paste0("run_kofamscan_parallel/",key,".kofamscan.tsv")
key.out <- paste0("reformatted_kegg_annotations-temp/",key,".kofamscan_reformatted.tsv")

jobfile <- paste("Rscript format_genomes_for_keggdecoder.R", key.in, key.out)

write.table(x = jobfile, file = "server-scripts/generated/2023-09-19_run_Keggdecoder/reformat_kofamscan_output.jobfile", row.names = F, col.names = F, quote = F)
