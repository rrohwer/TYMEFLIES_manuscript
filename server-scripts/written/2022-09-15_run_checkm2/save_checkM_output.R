# RRR
# I ran checkM2 per folder, so need to move all the folders off of scratch
# Also need to zip the protein_files folders before trying to transfer them

bins <- readRDS("data/2022-09-14_list_of_all_bins/bin_key.rds")

sample.names <- unique(bins$sample)


rename.script <- "#!/bin/bash"

for (s in sample.names){
  rename.script <- c(rename.script, "",
                     paste0("cd ",s),
                     paste0("mv protein_files ", s),
                     paste0("mv diamond_output/DIAMOND_RESULTS.tsv diamond_output/",s,".tsv"),
                     paste0("mv quality_report.tsv ",s,".tsv"),
                     "cd ../")
}


zip.script <- "#!/bin/bash"

for (s in sample.names){
  zip.script <- c(zip.script, "",
                  paste0("cd ",s),
                  paste0("tar cvfz ",s,".tar.gz ",s),
                  "cd ../")
}


org.script <- c("#!/bin/bash",
                "mkdir checkm2_protein_files",
                "mkdir checkm2_diamond_files",
                "mkdir checkm2_quality_files")

for (s in sample.names){
  org.script <- c(org.script, "",
                  paste0("cd ",s),
                  paste0("mv ",s,".tar.gz ../checkm2_protein_files"))
  
  org.script <- c(org.script,
                  paste0("mv diamond_output/",s,".tsv ../checkm2_diamond_files"))
  
  org.script <- c(org.script,
                  paste0("mv ",s,".tsv ../checkm2_quality_files"),
                  "cd ../")
}


write.table(x = rename.script, file = "server-scripts/generated/2022-09-15_run_checkm2/1-rename_checkm2_output.sh", append = F, quote = F, sep = "\n", row.names = F, col.names = F)
write.table(x = zip.script, file = "server-scripts/generated/2022-09-15_run_checkm2/2-zip_checkm2_output.sh", append = F, quote = F, sep = "\n", row.names = F, col.names = F)
write.table(x = org.script, file = "server-scripts/generated/2022-09-15_run_checkm2/3-organize_checkm2_output.sh", append = F, quote = F, sep = "\n", row.names = F, col.names = F)
