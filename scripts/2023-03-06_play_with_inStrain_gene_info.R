# RRR
# look at the gene.info table
# try to use data.table syntax, these datasets are getting big

library(data.table)
library(ggplot2)
genes <- fread(file = "~/Desktop/TYMEFLIES_large_data/instrain_example_files/ME2001-03-13pf_3300033984.IS/output/ME2001-03-13pf_3300033984.IS_gene_info.tsv", sep = "\t")

colnames(genes)

# what percent of genes are full vs partial?
x <- genes[ , .N, by = partial]
x$N[1] / (x$N[1] + x$N[2]) # 90% are full

# add a genome column

genes[ , genome := .(sub(pattern = "_scaffold_.*$",replacement = "", x = scaffold))]
colnames(genes)
head(genes)

# look at snvs/gene in different genomes, only full genes

my.genome <- unique(genes$genome)[2]
my.genome <- genes[genome == my.genome]

plot.data <- my.genome[partial == FALSE]
ggplot(data = plot.data, aes(y = genome, x = SNS_N_count))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_jitter(alpha = .2)

ggplot(data = genes, aes(x = SNV_N_count, y = nucl_diversity))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_point(alpha = .2)+
  scale_y_continuous(limits = c(0,.1))#+
  # scale_x_continuous(limits = c(0,300))
