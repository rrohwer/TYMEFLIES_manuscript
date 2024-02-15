# RRR

# ---- libraries and user input ----

library(data.table)

userprefs <- commandArgs(trailingOnly = TRUE)

genes <- userprefs[1] # ME2013-02-02s1D0_3300042325_group4_bin39_gene_info.tsv
cogdefs <- userprefs[2] # COG_descriptions.csv
anns <- userprefs[3] # ME2013-02-02s1D0_3300042325_group4_bin39.emapper.annotations
genomes <- userprefs[4] # drep_results_all_genomes_0.96.rds
threads <- userprefs[5]
data.folder <- userprefs[6]

# # COMMENT OUT TESTING PATHS
# cat("\n\n you left the mac paths defined! \n\n")
# genes <- "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39_gene_info.tsv"
# cogdefs <- "data/2023-04-30_COG_info/COG_descriptions.csv"
# anns <- "data/2023-04-19_test_merging_gene_info_with_annotations/ME2013-02-02s1D0_3300042325_group4_bin39.emapper.annotations"
# genomes <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"
# threads <- "1"
# data.folder <- "~/Desktop/test_plots/gene_info_data"
# summary.data.folder <- "~/Desktop/test_plots/selected_COGs_data"


# ---- read files and format folders ----

threads <- as.numeric(threads)
genes <- fread(file = genes, nThread = threads)
cogdefs <- fread(file = cogdefs, nThread = threads)
anns <- fread(file = anns, sep = "\t", quote = "", fill = T, skip = 4, nThread = threads)
genomes <- readRDS(genomes)

if(!dir.exists(file.path(data.folder))){
  dir.create(file.path(data.folder))
}

# ---- functions ----

get.fisher.p.value <- function(X){
  X <- matrix(X, nrow = 2, byrow = T)
  if(any(is.na(X))){
    X <- NA
  }else{
    X <- fisher.test(X)$p.value 
  }
  return(X)
}

add.fisher.p.value <- function(my.genes){
  my.mat <- as.matrix(my.genes[ ,.(SNV_N_count, SNV_S_count, SNS_N_count, SNS_S_count)])
  my.genes$MK.p.val <- apply(X = my.mat, MARGIN = 1, FUN = get.fisher.p.value)
  return(my.genes)
}

# ---- merge data ----

colnames(anns)[1] <- "query"
anns <- anns[-((nrow(anns) - 2):nrow(anns))] # stupid stats listed at end of file
anns <- anns[COG_category != "-", .(query, COG_category, Description)]

genes <- merge(x = genes, y = anns, by.x = "gene", by.y = "query", all.x = T, all.y = F) # all.y = F, can't end up with empty genome name if gene is missing in genes file 

genes <- merge(x = genes, y = cogdefs, by = "COG_category", all.x = T, all.y = F)

genes <- merge(x = genes, y = genomes, by.x = "genome", by.y = "bin.full.name", all.x = T, all.y = F)


# ---- find selected COGs ----

genes$mcdonald.kreitman <- genes$pNpS_variants / genes$dNdS_substitutions # < 1 positive selection, >1 negative or balancing selection

genes <- add.fisher.p.value(my.genes = genes)

# ---- save updated genes files ----

my.genome <- unique(genes$genome)

fwrite(x = genes, file = file.path(data.folder,paste0(my.genome,"_gene_info_combined_MK.tsv.gz")), sep = "\t", quote = F, row.names = F, append = F, nThread = threads)

cat("Finished",my.genome,"\n")


