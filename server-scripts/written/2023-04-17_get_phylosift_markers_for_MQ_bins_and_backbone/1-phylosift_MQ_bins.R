# RRR

# run phylosift on all MQ bins (not just Actinomycetia)
# cat output into alignments for each 96ANI group
# I don't think I need to worry about 30,000 files/folder because this is not on TACC

library(data.table)

all.bins <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
all.bins <- as.data.table(all.bins)

chosen.bins <- all.bins[!is.na(drep.cluster.0.96ANI), ]

# chosen.bins <- chosen.bins[1:5] # TEST

phylosift.jobs <- paste0("/home/baker/phylosift_v1.0.1/phylosift search --isolate --besthit /home/yggshare/current_projects/TYMEFLIES/tymeflies/metabat2_bins_renamed/", chosen.bins$tymeflies.name,"/",chosen.bins$bin.filename," && ",
                      "/home/baker/phylosift_v1.0.1/phylosift align --isolate --besthit /home/yggshare/current_projects/TYMEFLIES/tymeflies/metabat2_bins_renamed/", chosen.bins$tymeflies.name,"/",chosen.bins$bin.filename," && ",
                      "sed 's/\\.fna.*$//' PS_temp/",chosen.bins$bin.filename,"/alignDir/concat.updated.1.fasta > MQ_bin_phylosift_markers/",chosen.bins$bin.full.name,".faa") 

cat.jobs <- "#!/bin/bash"
for (b in unique(chosen.bins$drep.cluster.0.96ANI)){
  my.bins <- chosen.bins[drep.cluster.0.96ANI == b]
  my.rep <- my.bins[winner == TRUE, bin.full.name]
  my.bins <- my.bins[ ,bin.full.name]
  cat.jobs <- c(cat.jobs, 
                paste0("cat ",paste(paste0("MQ_bin_phylosift_markers/", my.bins,".faa"), collapse = " "), " > 96ANI_group_alignments/", my.rep,".faa"))
}         

write.table(x = phylosift.jobs, file = "server-scripts/generated/2023-04-17_get_phylosift_markers_for_MQ_bins_and_backbone/phylosift_MQ_bins.jobfile", col.names = F, row.names = F, quote = F)
write.table(x = cat.jobs, file = "server-scripts/generated/2023-04-17_get_phylosift_markers_for_MQ_bins_and_backbone/combine_96ANI_group_alignments.sh", col.names = F, row.names = F, quote = F)



