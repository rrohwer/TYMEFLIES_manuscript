# RRR

# for each sample folder, cat all the bin scaffolds into a single file in parallel
# then at the end cat all the sample-folder scaffolds into the same file

bin.key <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

sample.folders <- unique(bin.key$tymeflies.name)
  
grep.jobfile <- paste0("cat /home/yggshare/current_projects/TYMEFLIES/tymeflies/metabat2_bins_renamed/",sample.folders,"/* | grep \">\" > ",sample.folders,"_scaffolds.txt")

write.table(x = grep.jobfile, file = "server-scripts/generated/2023-04-12_get_contig_name_key/3-get_bin_scaffold_names.jobfile", quote = F, row.names = F, col.names = F)

# run as
# cat 3-get_bin_scaffold_names.jobfile | parallel -j 60