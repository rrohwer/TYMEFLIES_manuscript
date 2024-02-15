# RRR
# Run phylosift on the "backbone" genomes

library(data.table)

chosen.bins <- fread("data/2023-04-16_GTDB_Actinomycetia_backbone/5b_Actinomycetia_Backbone_Key.tsv")
bin.names <- sub(pattern = "\\.fna", replacement = "", x = chosen.bins$Name)

# chosen.bins <- chosen.bins[1:5] # TEST

phylosift.jobs <- paste0("/home/baker/phylosift_v1.0.1/phylosift search --isolate --besthit /home/yggshare/current_projects/TYMEFLIES/tymeflies/actinomycetia_tree/actinomycetia_backbone_genomes/", chosen.bins$Name," && ",
                         "/home/baker/phylosift_v1.0.1/phylosift align --isolate --besthit /home/yggshare/current_projects/TYMEFLIES/tymeflies/actinomycetia_tree/actinomycetia_backbone_genomes/", chosen.bins$Name," && ",
                         "sed 's/\\.fna.*$//' PS_temp/",chosen.bins$Name,"/alignDir/concat.updated.1.fasta > Actinomycetia_backbone_markers/",bin.names,".faa") 

     

write.table(x = phylosift.jobs, file = "server-scripts/generated/2023-04-17_get_phylosift_markers_for_MQ_bins_and_backbone/phylosift_Actinomycetia_backbone.jobfile", col.names = F, row.names = F, quote = F)
