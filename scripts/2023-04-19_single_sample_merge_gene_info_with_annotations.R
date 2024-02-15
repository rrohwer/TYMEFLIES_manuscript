# RRR
# Try for just one bin

library(data.table)

gene <- fread("data/2023-04-19_test_merging_gene_info_with_annotations/ME2000-03-15pf_3300044539_group1_bin132_gene_info.tsv")
gff <-  fread("data/2023-04-19_test_merging_gene_info_with_annotations/ME2000-03-15pf_3300044539_group1_bin132.gff3")
tsv <-  fread("data/2023-04-19_test_merging_gene_info_with_annotations/ME2000-03-15pf_3300044539_group1_bin132.tsv", fill = T)

gff.comments <- fread("data/2023-04-19_test_merging_gene_info_with_annotations/comment.lines.txt", skip = 3)

boxplot(x = list("instrain" = gene$end - gene$start, "interproscan" = gff.comments$V4 - gff.comments$V3)) # aha! 3X as much???
summary(gff.comments$V4)


# gff is base-1, instrain is base-0
# gff.comments[ ,V3 := V3 -1]
# gff.comments[ ,V4 := V4 -1]

# interproscan counted amino acids, instrain counted nucleotides
gff.comments[ ,V4 := V4 * 3]
# gff.comments[ ,V3 := V3 * 3] # DON'T multiply the START by 3, only the END... but some starts should be x3 if they're not 1?

boxplot(x = list("instrain" = gene$end - gene$start, "interproscan" = gff.comments$V4 - gff.comments$V3)) # aha! 3X as much???

# AHHHHHHHHHH but still not very many
x <- merge(x = gene, y = gff.comments, by.x = c("gene","start","end"), by.y = c("V2","V3","V4"))

# most (but not all) of the genes exist in both datasets
x <- merge(x = gff.comments, y = gene, by.x = c("V2"), by.y = c("gene"))


y <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
sample.names <- unique(y$tymeflies.name)

my.script <- paste0("cat ", sample.names, " > ", "new")

write.table( quote = F, col.names = F, row.names = F)


# instrain take in prodigal genes.fna file
# interproscan takes in prodigal genes.faa file
# but it could read in the fna file instead- but then it calls the ORF itself, won't be using prodigal right?

# what about the stars I deleted? well every gene didn't have a star, so some should have still matched
