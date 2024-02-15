# RRR
# when filter to only significant KOs, there are basically no complete pathways in the genome!
# need to choose how to decide the cutoffs better.

# ---- set-up ----

library(data.table)

pathways.key <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/ko_pathway_links.txt", header = F)

pathways.names <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/pathway_list.txt", header = F)

ko.names <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/ko_list.txt", header = F, sep = "\t")

genome.annotations <- fread(file = "data/2023-09-13_KEGG_annotations/ME2011-09-21_3300043464_group3_bin69.kofamscan.tsv")

# ---- format files ----

colnames(pathways.key) <- c("ko","path")
pathways.key[ ,`:=`(ko = sub(pattern = "ko:", replacement = "", x = ko),
                      path = sub(pattern = "path:", replacement = "", x = path))]
pathways.key <- pathways.key[grep(pattern = "map", x = path)]

colnames(pathways.names) <- c("path","pathway.description")

colnames(ko.names) <- c("ko","ko.description")

genome.annotations <- genome.annotations[-1]
colnames(genome.annotations) <- c("signif","gene","ko","threshold","score","e.value","ko.description")


# ---- which pathways are in the genome ----

# for each pathway, list its KOs 

pathways.list <- list(NULL)
for (p in unique(pathways.key$path)){
  k <- pathways.key[path == p]
  k <- list(k$ko)
  names(k) <- p
  pathways.list <- c(pathways.list, k)
}
pathways.list

# for each pathway, how many of its KOs are present?

genome.annotations <- genome.annotations[signif == "*"]

genome.pathways <- pathways.key[ , .N, by = .(path)]
genome.pathways$genome.N <- 0

for (p in genome.pathways$path){
  path.kos <- pathways.list[[p]]
  x <- intersect(x = path.kos, genome.annotations$ko)
  genome.pathways[path == p, genome.N := length(x)]
}

genome.pathways[ ,perc.complete := genome.N / N * 100]

# put names to the pathways and sort by how big they are

genome.pathways <- merge(x = pathways.names, y = genome.pathways, by = "path", all.x = F, all.y = T)

genome.pathways <- genome.pathways[order(N, decreasing = T)]


# ---- which pathways do the genome's KOs belong to ----

# only include genes with known annotation and known pathways

genome.annotations <- merge(x = genome.annotations, y = pathways.key, by = "ko", all.x = F, all.y = F)

# only include pathways that are at least 50% complete

genome.annotations <- merge(x = genome.annotations, y = genome.pathways[perc.complete > 10], by = "path", all.x = F, all.y = F)


# ---- pull in genes under selection data ----

library(lubridate)
library(ggplot2)

genes <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz")
genes <- genes[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman)]

genes <- merge(x = genome.annotations, y = genes, by = "gene", all.x = F, all.y = T, allow.cartesian = T)

colnames(genes)[colnames(genes) == "N"] <- "N.pathway"
genes <- genes[ , .N, by = .(date, pathway.description, N.pathway)]
genes <- genes[order(N.pathway, decreasing = T)]

# genes <- dcast(data = genes, formula = pathway.description ~ date, value.var = "N")
# genes <- as.matrix(genes, rownames = "pathway.description")

# my.colors <- colorRampPalette(colors = c("white","red4"))

genes[ ,date := parse_date_time(date, orders = "ymd")]

pdf(file = "~/Desktop/test69_all_paths.pdf", width = 24, height = 24)
ggplot(data = genes, aes(x = date, y = pathway.description, color = N))+
  geom_point(aes(size = N), alpha = .5)+
  theme_bw()+
  scale_color_gradient(low = "grey90", high = "red4")
dev.off()

# try removing the huge categories...
genes <- genes[pathway.description != "Metabolic pathways" & 
                 pathway.description != "ABC transporters" & 
                 pathway.description != "Two-component system" &
                 pathway.description != "Microbial metabolism in diverse environments" &
                 pathway.description != "Biosynthesis of secondary metabolites"]

pdf(file = "~/Desktop/test69_some_paths.pdf", width = 24, height = 24)
ggplot(data = genes, aes(x = date, y = pathway.description, color = N))+
  geom_point(aes(size = N), alpha = .5)+
  theme_bw()+
  scale_color_gradient(low = "grey90", high = "red4")
dev.off()
