# RRR
# read in the per-genome SNV files
# remove dates where the genome is at too low abundance for reliable SNV calling
# create a distance matrix based on the percent reference allele
# create an NMDS object from the distance table
# do that for both ALL snvs and only the nonsynonymous SNVs

# NOTE: create the 4 output folders before running the script
# NOTE: the /tmp folder filled up and caused an error with 118 threads, deleted the Rtemp files and ran with 50 and it didn't throw the error.

# example syntax:
#  Rscript 1-calculate_distance_matrix_and_NMDS.R per-genome_SNVs/ME2011-09-21_3300043464_group3_bin69_SNVs.tsv.gz multivariate_SNV_analysis/distance_matrices_all_SNVs multivariate_SNV_analysis/distance_matrices_nonsynon_SNVs multivariate_SNV_analysis/NMDS_objects_all_SNVs multivariate_SNV_analysis/NMDS_objects_nonsynon_SNVs 1

# ---- set up ----

library(data.table)
library(vegan)

userprefs <- commandArgs(trailingOnly = TRUE)
per.genome.snv.file <- userprefs[1]
genome.info.file <- userprefs[2]
output.dist.folder.all <- userprefs[3]
output.dist.folder.N <- userprefs[4]
output.nmds.object.folder.all <- userprefs[5]
output.nmds.object.folder.N <- userprefs[6]
num.threads <- as.numeric(userprefs[7])

# # local path testing
# per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-21_3300043464_group3_bin69_SNVs.tsv.gz"
# genome.info.file <- "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz"
# output.dist.folder.all <- "data/2023-11-01_multidimensional_SNV_analysis"
# output.dist.folder.N <- "data/2023-11-01_multidimensional_SNV_analysis"
# output.nmds.object.folder.all <- "data/2023-11-01_multidimensional_SNV_analysis"
# output.nmds.object.folder.N <- "data/2023-11-01_multidimensional_SNV_analysis"
# num.threads <- 1

# ---- import data ----

my.genome <- sub(pattern = "^.*/", replacement = "", x = per.genome.snv.file)
my.genome <- sub(pattern = "_SNVs.*$", replacement = "", x = my.genome)
cat("Processing",my.genome,"\n")

genome.info <- fread(file = genome.info.file, nThread = num.threads)
genome.info <- genome.info[genome == my.genome]

x <- "character"
names(x) <- "date"
snv.instrain <- fread(file = per.genome.snv.file, nThread = num.threads, colClasses = x)

# ---- process data ----

# Remove samples with coverage too low to reliably call SNVs
genome.info <- genome.info[ ,above.breadth.cutoff := (breadth / breadth_expected) >= .5]
genome.info <- genome.info[above.breadth.cutoff == TRUE & coverage_median > 10, sample]

if (length(genome.info) > 2){
  snv.instrain <- snv.instrain[sample %in% genome.info]
  
  # When the SNV doesn't exist in one sample, call it 1 because if the SNV is not there then 100% of the reads are the reference
  snvs <- dcast(data = snv.instrain, formula = scaffold + position + gene + mutation_type ~ sample, value.var = "ref_freq", fill = 1)
  
  snvs[ ,snv.num := paste0("snv",1:nrow(snvs))] # need unique row names
  
  # one version of data with all SNVs: N, S, and not even on a gene
  snv.mat.all <- as.matrix(snvs[ ,-c(1:4)], rownames = "snv.num")
  snv.mat.all <- t(snv.mat.all)
  
  # one verion with only the nonsynonymous snvs 
  snv.mat.N <- as.matrix(snvs[mutation_type == "N", -c(1:4)], rownames = "snv.num")
  snv.mat.N <- t(snv.mat.N)
  
  # use euclidean distance, because unlike species with an ecosystem carrying capacity (so you use bray curtis) the SNVs do not have a constant max
  snv.dist.euc.all <- vegdist(x = snv.mat.all, method = "euclidean")
  snv.nmds.euc.all <- metaMDS(comm = snv.dist.euc.all, trymax = 1000)
  
  snv.dist.euc.N <- vegdist(x = snv.mat.N, method = "euclidean")
  snv.nmds.euc.N <- metaMDS(comm = snv.dist.euc.N, trymax = 1000)
  
  # ---- export data ----
  
  # save as rds objects because then they will read in correctly as matrices and lists
  
  saveRDS(object = snv.dist.euc.all, file = file.path(output.dist.folder.all, paste0(my.genome,"_all_SNV_euclidean_distance_matrix.rds")))
  
  saveRDS(object = snv.dist.euc.N, file = file.path(output.dist.folder.N, paste0(my.genome,"_nonsynonymous_SNV_euclidean_distance_matrix.rds")))
  
  saveRDS(object = snv.nmds.euc.all, file = file.path(output.nmds.object.folder.all, paste0(my.genome,"_all_SNV_euclidean_distance_NMDS_object.rds")))
  
  saveRDS(object = snv.nmds.euc.N, file = file.path(output.nmds.object.folder.N, paste0(my.genome,"_nonsynonymous_SNV_euclidean_distance_NMDS_object.rds")))
  
}

