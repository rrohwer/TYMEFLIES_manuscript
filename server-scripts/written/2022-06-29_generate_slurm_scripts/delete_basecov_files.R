# RRR
# for first group, generated the bbmap basecov stats
# but it is too big to save, so delete them all with this script

x <- readRDS(file = "data/2022-06-27_binning_groups/mapping_groups_list.rds")
output.folder <- "server-scripts/generated/2022-06-29_slurm_scripts/"

group <- names(x)[1]
for (group in names(x)[1:2]){ # just for test and group1
  
  sample.names <- x[[group]]
  
  delete.script <- "#!/bin/bash"
  
  for (assembly in sample.names){
    delete.script <- c(delete.script, 
                       paste0("cd ", assembly),
                       "rm *.basecov",
                       "cd ../")
  }
  
  write.table(x = delete.script, file = file.path(output.folder, group, "delete_basecov.sh"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
}
