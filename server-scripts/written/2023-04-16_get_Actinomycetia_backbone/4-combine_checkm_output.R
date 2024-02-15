# so the way I ran it was kind of annoying, made a separate folder for each bin....

bin.names <- read.table("data/2023-04-16_GTDB_accessions/3a_gtdb_downloaded_bin_names.txt")
bin.names <- bin.names$V1
folder.names <- sub(pattern = "\\.fna\\.gz$", replacement = "", x = bin.names)

script <- c("#!/bin/bash",
            paste0("cat output/",folder.names[1],"/quality_report.tsv > checkm2_output_GTDB_Actinomycetia.tsv"),
            paste0("cat output/",folder.names[-1],"/quality_report.tsv >> checkm2_output_GTDB_Actinomycetia.tsv"))

write.table(x = script, file = "server-scripts/generated/2023-04-16_get_Actinomycetia_backbone/4-combine_checkm2_output.sh", quote = F, row.names = F, col.names = F)