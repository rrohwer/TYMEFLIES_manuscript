# RRR
# have a folder called dRep_ANIrange_output
# for each dRep ANI
#   cp the wanted output to the new folder, changing it's name to include the ANI

ani.range <- seq(from = .90, to = .99, by = .01)


copy.genome.information <- paste0("cp drep_output_",ani.range,"/data_tables/genomeInformation.csv dRep_ANIrange_output/genomeInformation_",ani.range,".csv")

copy.cdb <- paste0("cp drep_output_",ani.range,"/data_tables/Cdb.csv dRep_ANIrange_output/Cdb_",ani.range,".csv")

copy.wdb <- paste0("cp drep_output_",ani.range,"/data_tables/Wdb.csv dRep_ANIrange_output/Wdb_",ani.range,".csv")


shell.script <- c("#!/bin/bash",
                  copy.genome.information,
                  copy.cdb,
                  copy.wdb)

write.table(x = shell.script, file = "server-scripts/generated/2023-02-24_run_dRep_range_of_ANIs/copy_dRep_output_files.sh",
            row.names = F, col.names = F, quote = F)
