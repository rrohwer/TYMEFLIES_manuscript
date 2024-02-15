# RRR

bin.key <- readRDS("data/2022-09-14_list_of_all_bins/bin_key.rds")

shell.script <- "#!/bin/bash"

shell.script <- c(shell.script,
                  paste0("cd ", bin.key$sample, " ; ",
                       "mv ", bin.key$bin, " ", bin.key$name, " ; ",
                       "cd ../"))
shell.script[1:5]

write.table(x = shell.script, file = "server-scripts/generated/2022-09-15_rename_bins/rename_bins.sh", row.names = F, quote = F, col.names = F, sep = "\n")
