# RRR

library(data.table)
library(lubridate)
library(ggplot2)
library(ggtree)

genes <- fread("data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.consistently_selected_genes-by_sample.tsv.gz")

# genes <- genes[ ,.(gene,ko,ko.description,path,pathway.description, module, module.description, sample,date,year,yday,season,genome.pres, mcdonald.kreitman,MK.p.val,pres,pos)]

# ---- format data ----

# format table
genes[ ,date := parse_date_time(x = date, orders = "ymd")]
genes[ ,pathway.description := gsub(pattern = "; ", replacement = "\n", x = pathway.description)]
genes[ ,module.description := gsub(pattern = "; ", replacement = "\n", x = module.description)]

# only call it pos selection if p-value is <= .25, otherwise p-value is NA
genes$positive.selection.pvalue <- as.numeric(1)
genes[mcdonald.kreitman < 1 & MK.p.val <= .25, positive.selection.pvalue := MK.p.val]

# subset to overall
genes.group <- genes[consistent.in == "overall"]

# sort based on occurrence patterns over time
genes.wide <- dcast(data = genes.group, formula = gene + ko + ko.description + module + module.description + path + pathway.description ~ sample, value.var = "positive.selection.pvalue")
genes.mat <- as.matrix(x = genes.wide[ ,-c(2:7)], rownames = 1)
genes.dist <- dist(x = genes.mat, method = "euclidean", diag = T)
genes.clust <- hclust(d = genes.dist, method = "ward.D")
genes.dend <- as.dendrogram(genes.clust)


# ---- make plots ----

# ----
pdf(file = "figures/2023-09-25_exploring_KOs_and_pathways/clustered_heatmap-By_gene.pdf", width = 24, height = 8)
p.tree <- ggtree(genes.dend) %<+% 
  genes.wide +
  geom_label(aes(label = ko.description), label.size = NA, hjust = 0)
gheatmap(p = p.tree, data = genes.mat, offset = 3, low = "red4", high = adjustcolor("grey70",.25), colnames = F)+
  theme(plot.margin = margin(0,0,0,0))+
  labs(title = "Gene KO annotations")
dev.off()

# ----
# ----
pdf(file = "figures/2023-09-25_exploring_KOs_and_pathways/clustered_heatmap-By_pathway.pdf", width = 24, height = 26)
p.tree <- ggtree(genes.dend) %<+% 
  genes.wide +
  geom_label(aes(label = pathway.description), label.size = NA, hjust = 0)
gheatmap(p = p.tree, data = genes.mat, offset = 2, low = "red4", high = adjustcolor("grey70",.3), colnames = F)+
  theme(plot.margin = margin(0,0,0,0))+
  labs(title = "Gene pathway options")
dev.off()

# ----

# ----
pdf(file = "figures/2023-09-25_exploring_KOs_and_pathways/clustered_heatmap-By_module.pdf", width = 24, height = 8)
p.tree <- ggtree(genes.dend) %<+% 
  genes.wide +
  geom_label(aes(label = module.description), label.size = NA, hjust = 0)
gheatmap(p = p.tree, data = genes.mat, offset = 2.5, low = "red4", high = adjustcolor("grey70",.3), colnames = F)+
  theme(plot.margin = margin(0,0,0,0))+
  labs(title = "Gene module options")
dev.off()

# ----


# and for some reason this doesn't work with gene as the column name WTF
genes.wide$new.col <- sub(pattern = "ME2011-09-21_3300043464_group3_bin69_", replacement = "", x = genes.wide$gene)
# ----
pdf(file = "figures/2023-09-25_exploring_KOs_and_pathways/clustered_heatmap-By_gene_name.pdf", width = 24, height = 8)
p.tree <- ggtree(genes.dend) %<+% 
  genes.wide +
  geom_label(aes(label = new.col), label.size = NA, hjust = 0)
gheatmap(p = p.tree, data = genes.mat, offset = 1, low = "red4", high = adjustcolor("grey70",.3), colnames = F)+
  theme(plot.margin = margin(0,0,0,0))+
  labs(title = "Gene names")
dev.off()
# ----


genes.wide[ ,1:3]
genes.wide[grep(pattern = "pyrimidine", x = pathway.description, ignore.case = T),1:3]





