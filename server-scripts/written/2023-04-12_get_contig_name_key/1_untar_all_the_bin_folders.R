# RRR

# need to untar all the bin folder on the Naz
# make a simple job file and just run it in parallel

folder.names <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

folder.names <- unique(folder.names$tymeflies.name)

unzip.jobfile <- paste0("tar xvf ",folder.names,".tar")

write.table(x = unzip.jobfile, file = "server-scripts/generated/2023-04-12_get_contig_name_key/1_untar_all_the_bin_folders.jobfile", quote = F, row.names = F, col.names = F)

# run as
# cat 1_untar_all_the_bin_folders.jobfile | parallel -j 60

# ----

# also need to change the file permissions to read-only because I don't want to edit the bin files by accident

chmod.jobfile <- paste0("find ",folder.names,"/. -exec chmod -w {} ","\\\\",";")
cat(chmod.jobfile)

write.table(x = chmod.jobfile, file = "server-scripts/generated/2023-04-12_get_contig_name_key/2_change_the_fucking_file_permissions.jobfile", quote = F, row.names = F, col.names = F)

# and then find \\ and replace with \ 
# run as
# cat 2_change_the_fucking_file_permissions.jobfile | parallel -j 60
