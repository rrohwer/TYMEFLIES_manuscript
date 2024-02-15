# RRR

# combine the heatmap consistent gene keys into a large file
# this will let me compare them
# not having the per-sample info should help it not be too big
# think I can do on my computer only

# no wait, need to do on the server, will be more convenient b/c can put all selected genes in 1 table

# example server call (took less than a minute)
# Rscript combine_consistent_gene_tables.R consistent_genes_analysis/key_genes kegg_annotations_with_pathways consistent_genes_analysis/key_genes_all_genomes-anns_separate_lines.tsv.gz 45 

# ---- set-up ----

library(data.table)

userprefs <- commandArgs(trailingOnly = TRUE)
folder.path <- userprefs[1]
folder.path.2 <- userprefs[2]
output.file <- userprefs[3]
threads <- as.numeric(userprefs[4])

# # example local paths
# folder.path <- "data/2023-10-16_consistent_genes_example_data/consistent_genes"
# folder.path.2 <- "data/2023-10-16_consistent_genes_example_data/gene_annotations/"
# output.file <- "data/2023-10-16_consistent_genes_example_data/combined_file/consistently_selected_genes_KOs_separate_example_genomes.tsv.gz"
# threads <- 1

# ---- combine files ----

my.files <- list.files(folder.path, full.names = F, pattern = "^ME")
my.files <- data.table("consistent.file" = my.files, "genome" = sub(pattern = ".consistently_selected_genes.tsv.gz", replacement = "", x = my.files))

my.files.2 <- list.files(folder.path.2, full.names = F, pattern = "^ME")
my.files.2 <- data.table("expanded.file" = my.files.2, "genome" = sub(pattern = ".ko-pathways_anns-separate.tsv.gz", replacement = "", x = my.files.2))

my.files <- merge(x = my.files, y = my.files.2, by = "genome", all = T)

table.list <- list()
for (r in 1:nrow(my.files)){
  genome.name <- my.files[r, genome]
  cat("reading genome",genome.name,"\n")
  consistent <- fread(file = file.path(folder.path, my.files[r, consistent.file]), nThread = threads)
  expanded <- fread(file = file.path(folder.path.2, my.files[r, expanded.file]), nThread = threads)
  temp <- merge(x = consistent[ ,c("gene","consistent.in","pres","pos","pos.perc","is.Q4","is.outlier")],
                y = expanded[ ,c("gene","threshold","bit.score","e.value","signif","relaxed.signif.6","ko","ko.description","module","module.description","module.N","genome.module.N","module.perc.complete","path","pathway.description","pathway.N","genome.pathway.N","pathway.perc.complete")],
                by = "gene", all = T)
  temp$genome <- genome.name
  setcolorder(x = temp, neworder = c("genome","gene","consistent.in","pres","pos","pos.perc","is.Q4","is.outlier",
                                     "threshold","bit.score","e.value","signif","relaxed.signif.6",
                                     "ko","ko.description","module","module.description","module.N","genome.module.N","module.perc.complete","path","pathway.description","pathway.N","genome.pathway.N","pathway.perc.complete"))
  table.list <- c(table.list, list(temp))
}
big.table <- rbindlist(table.list)
colnames(big.table)

cat("\nwriting combined file",output.file,"\n")
fwrite(x = big.table, file = output.file, nThread = threads)

