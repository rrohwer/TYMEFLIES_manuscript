# RRR
# make a jobfile to run with parallel
# cat get_SNV_tables.jobfile | parallel -j 100 

key <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
key <- key[key$winner, ]

file.list <- paste0("per-genome_SNVs/", key$bin.full.name,"_SNVs.tsv.gz")

jobfile <- paste0("Rscript get_SNV_pres_abs_tables.R ", file.list, " 1 SNV_pres_abs_tables new_SNV_stats_tables")

write.table(x = jobfile, file = "server-scripts/generated/2023-07-18_get_SNV_tables/get_SNV_tables.jobfile", quote = F, row.names = F, col.names = F)

# nevermind brett's servers are fucked do it on tacc

launcherfile <- paste0("module load Rstats/4.0.3 ; Rscript get_SNV_pres_abs_tables.R ", file.list, " 1 SNV_pres_abs_tables new_SNV_stats_tables")

write.table(x = launcherfile, file = "server-scripts/generated/2023-07-18_get_SNV_tables/get_SNV_tables.launcher", quote = F, row.names = F, col.names = F)

# fuuuuuuck tacc has some accounting error. just run it for 27 hrs instead of 16 min fuck this

shellscript <- c("#!/bin/bash", jobfile)

write.table(x = shellscript, file = "server-scripts/generated/2023-07-18_get_SNV_tables/get_SNV_tables.sh", quote = F, row.names = F, col.names = F)
