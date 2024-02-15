# RRR
# use parallel to run in parallel on midgard
# make a jobs file very similar to a launcher file
# run it with 
# cat jobsfile | parallel -j numjobs

all.bins <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
chosen.bins <- all.bins[all.bins$winner & all.bins$phylum == "p__Actinobacteriota" & !is.na(all.bins$phylum), ]

jobs.file <- c(paste0("/home/baker/phylosift_v1.0.1/phylosift search --isolate --besthit ../../jyuan/drep_genomes_96ANI/",chosen.bins$bin.filename," && ",
                      "/home/baker/phylosift_v1.0.1/phylosift align --isolate --besthit ../../jyuan/drep_genomes_96ANI/",chosen.bins$bin.filename," && ",
                      "sed 's/\\.fna.*$//' PS_temp/",chosen.bins$bin.filename,"/alignDir/concat.updated.1.fasta > ",paste0(chosen.bins$bin.full.name,".faa")))



# cat(paste0("sed 's/\\.fna.*$//'")) # note that R requires an extra backslash to escape the backslash, but this writes correctly as sed 's/\.fna.*$//'

write.table(x = jobs.file, file = "server-scripts/generated/2023-03-13_run_phylosift_on_midgard/phylosift_actinos.jobs", quote = F, row.names = F, col.names = F)


