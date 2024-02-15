
sample.names <- readRDS("data/2022-06-27_binning_groups/2022-06-27_bigtyme.rds")

sample.names <- sample.names$sample.names

sample.names

tar.gz <- paste0(sample.names, ".tar.gz")
tar <- paste0(sample.names, ".tar")

script.tar.gz <- c("#!/bin/bash", paste("tar cvfz", tar.gz, sample.names))
script.tar <- c("#!/bin/bash", paste("tar cvf", tar.gz, sample.names))
script.untar.gz <- c("#!/bin/bash", paste("tar xvfz", tar.gz))
script.untar <- c("#!/bin/bash", paste("tar xvf", tar.gz))


write.table(x = script.tar.gz, file = "server-scripts/generated/2022-11-01_zip_sample_folders/tar_gz_all_sample_folders.sh", quote = F, row.names = F, col.names = F)
write.table(x = script.tar, file = "server-scripts/generated/2022-11-01_zip_sample_folders/tar_all_sample_folders.sh", quote = F, row.names = F, col.names = F)
write.table(x = script.untar.gz, file = "server-scripts/generated/2022-11-01_zip_sample_folders/untar_gz_all_sample_folders.sh", quote = F, row.names = F, col.names = F)
write.table(x = script.untar, file = "server-scripts/generated/2022-11-01_zip_sample_folders/untar_all_sample_folders.sh", quote = F, row.names = F, col.names = F)
