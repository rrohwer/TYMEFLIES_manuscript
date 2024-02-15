x <- readRDS(file = "data/2022-09-14_list_of_all_bins/bin_key.rds")
head(x)

sample.names <- unique(x$sample)

generated.script.folder <- "server-scripts/generated/2022-09-22_run_gtdb-tk/"

shell.script <- "#!/bin/bash"
for (s in sample.names){
  shell.script <- c(shell.script,
                    paste0("tar cvf ",s,".tar ",s))
  
}

write.table(x = shell.script, file = file.path(generated.script.folder,"tar_all_the_bins.sh"), quote = F, sep = "\n", row.names = F, col.names = F)


# OK also need to untar the bins

shell.script <- "#!/bin/bash"
for (s in sample.names){
  shell.script <- c(shell.script,
                    paste0("tar xvf ",s,".tar"))
  
}

write.table(x = shell.script, file = file.path(generated.script.folder,"untar_all_the_bins.sh"), quote = F, sep = "\n", row.names = F, col.names = F)

# or could just do 
# find . -exec tar xvf {} \;