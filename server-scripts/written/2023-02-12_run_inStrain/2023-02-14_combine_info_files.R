# RRR
# combine the output files into folders of just that type of output. 
# Then I can zip and transfer just the one file type
# this will make it easier to combine the files into one spreadsheet.


all.bins <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")

sample.names <- unique(all.bins$tymeflies.name)

genome.info.files <- paste0(sample.names, ".IS_genome_info.tsv")

cp.genome.files.sh <- c("#!/bin/bash",
                        paste0("cp ",sample.names,".IS/output/",genome.info.files," genome_info_files/"))

write.table(x = cp.genome.files.sh, file = "server-scripts/generated/2023-02-12_run_inStrain/copy_genome_info_files_to_folder.sh",
            quote = F, row.names = F, col.names = F)



gene.info.files <- paste0(sample.names, ".IS_gene_info.tsv")

cp.gene.files.sh <- c("#!/bin/bash",
                        paste0("cp ",sample.names,".IS/output/",gene.info.files," gene_info_files/"))

write.table(x = cp.gene.files.sh, file = "server-scripts/generated/2023-02-12_run_inStrain/copy_gene_info_files_to_folder.sh",
            quote = F, row.names = F, col.names = F)
