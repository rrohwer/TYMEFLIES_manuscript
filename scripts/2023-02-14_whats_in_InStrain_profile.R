# RRR
# OK, I have some inStrain results from inStrain profile command
# but the built-in plotting all failed due to python version being > 3.8
# Also, there was a timeout error, but it looks like it completed and maybe just didn't send a correct exit code?
# so anyway, what files do I have?

library(data.table)

gene.info <- "~/Desktop/ME2000-03-15pf_3300044539.IS_gene_info.tsv"
gene.info <- fread(file = gene.info, sep = "\t", header = T)
colnames(gene.info)
dim(gene.info)
gene.info
gene.info$genome <- sub(pattern = "_scaffold_.*$",replacement = "", x = gene.info$scaffold)
unique(gene.info$genome)

genome.info <- "~/Desktop/TYMEFLIES_large_data/instrain_example_files/ME2001-03-13pf_3300033984.IS/output/ME2001-03-13pf_3300033984.IS_genome_info.tsv"
genome.info <- read.table(file = genome.info, sep = "\t", header = T)
colnames(genome.info)
dim(genome.info)

# this is the 4GB file
linkage <- "~/Desktop/TYMEFLIES_large_data/instrain_example_files/ME2001-03-13pf_3300033984.IS/output/ME2001-03-13pf_3300033984.IS_linkage.tsv" 
linkage <- read.table(file = linkage, header = TRUE, sep = "\t") # too long, open in bbedit faster

mapping.info <- "~/Desktop/TYMEFLIES_large_data/instrain_example_files/ME2001-03-13pf_3300033984.IS/output/ME2001-03-13pf_3300033984.IS_mapping_info.tsv" 
mapping.info <- read.table(file = mapping.info, header = T, sep = "\t")
colnames(mapping.info)
dim(mapping.info)

scaffold.info <- "~/Desktop/TYMEFLIES_large_data/instrain_example_files/ME2001-03-13pf_3300033984.IS/output/ME2001-03-13pf_3300033984.IS_scaffold_info.tsv" 
scaffold.info <- read.table(file = scaffold.info, header = T, sep = "\t")
colnames(scaffold.info)

# this one is 1 GB
snvs <- "~/Desktop/TYMEFLIES_large_data/instrain_example_files/ME2001-03-13pf_3300033984.IS/output/ME2001-03-13pf_3300033984.IS_SNVs.tsv" 
snvs <- read.table(file = snvs, header = T, sep = "\t")
colnames(snvs)

