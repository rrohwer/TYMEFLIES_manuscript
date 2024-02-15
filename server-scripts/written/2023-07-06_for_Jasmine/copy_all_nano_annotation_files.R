# RRR
# 7.6.23 
# Jasmine needs a folder with just the Actinomycetia annotation files to make her plot. 
# Make a shell script to copy just those files

key <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

key <- key[key$winner, ]
key <- key[key$order == "o__Nanopelagicales" & !is.na(key$order), ]

# name format: ME2019-04-24_3300042353_group8_bin52.emapper.annotations.gz
ann.file <- paste0(key$bin.full.name,".emapper.annotations.gz")

shell.script <- c("#!/bin/bash",
                  paste0("cp ../../yggshare/current_projects/TYMEFLIES/tymeflies/eggnog_annotations/",ann.file," ./"))

write.table(x = shell.script, file = "server-scripts/generated/2023-07-06_for_Jasmine/copy_all_nano_annotation_files.sh", quote = F, row.names = F, col.names = F)
