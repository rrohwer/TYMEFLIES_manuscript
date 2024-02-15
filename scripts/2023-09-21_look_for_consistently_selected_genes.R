# RRR

library(data.table)
library(ggplot2)
library(lubridate)
library(patchwork)

genes <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz")

coverm <- fread(input = "data/2023-03-15_coverM_on_drep96/output/coverM_drep96_minANI93.txt.gz")
genome.info <- fread(input = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
my.genome <- "ME2011-09-21_3300043464_group3_bin69"


# list whether genome is there ----

colnames(coverm)[3] <- "rel.abund"
coverm <- coverm[ ,1:3]
coverm <- coverm[Genome == my.genome]
coverm[ ,date := substr(Sample, start = 3, stop = 12)]
coverm[ ,date := parse_date_time(x = date, orders = "ymd")]
coverm[ ,Sample := NULL]
coverm <- coverm[ ,.(rel.abund = mean(rel.abund)), by = date] # should choose unfiltered for the one dup day with ww and pf

genome.info <- genome.info[genome == my.genome]
genome.info[ ,date := parse_date_time(x = date, orders = "ymd")]
genome.info[ ,above.breadth.cutoff := (breadth / breadth_expected) >= .5, ]

genome <- merge(x = genome.info, y = coverm, by = "date")
genome <- genome[ ,.(genome, sample, rel.abund, above.breadth.cutoff)]
summary(genome$rel.abund[genome.info$above.breadth.cutoff])
genome$genome.pres <- FALSE
genome[rel.abund > .1 & above.breadth.cutoff == TRUE, genome.pres := TRUE]

# list whether gene is there and whether it's selected for ----

genes$pos <- FALSE
genes[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman), pos := TRUE]

genes$pres <- FALSE
summary(genes$coverage)
genes[coverage >= 5 & breadth >= 1, pres := TRUE]

genes[ , sample := sub(pattern = "\\.IS_gene_info\\.tsv", replacement = "", x = sample)]
genes[ , sample := sub(pattern = "\\.gz", replacement = "", x = sample)]

genes <- merge(x = genes, genome[ , .(sample, genome.pres)], by = "sample")

# how often are they selected vs. there?
consistency <- genes[ ,.(gene, sample, year, season, genome.pres, pres, pos)]

con.overall <- consistency[genome.pres == TRUE, .(pres = sum(pres), pos = sum(pos)), by = .(gene)]

con.pre <- consistency[year < 2012 & genome.pres == TRUE, .(pres = sum(pres), pos = sum(pos)), by = .(gene)]
con.post <- consistency[year >= 2012 & genome.pres == TRUE, .(pres = sum(pres), pos = sum(pos)), by = .(gene)]

con.pre.post <- merge(con.pre, con.post, by = "gene", suffixes = c(".pre",".post"))

# pull out outlier "key genes" specifically

key.genes <- con.overall[pos > 0]

outliers <- boxplot(key.genes[ ,pos]) # could do Q3 or something instead also to include more
key.genes <- key.genes[pos >= min(outliers$out)]

key.genes.info <- merge(x = key.genes, y = genes, by = "gene", all.x = TRUE, all.y = FALSE)

# ---- Look at stuff ----

# how many genes are inconsistently there? And are they the same ones coming under selection?
ggplot(data = con.overall, aes(y = pres))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(x = 0, color = pos), width = .36)+
  scale_color_gradient(low = adjustcolor("grey80",alpha = .2), high = "blue", name = "Times under\npositive selection")+
  ylab(label = "Times Gene is present\n(when genome is at least .1% abundance and 50% covered")+
  labs(title = "Are the genes under selection consistently in the genome?", 
       subtitle = "each dot is a gene, the color is how many times it was selected. all genes are included")


# Are there clear outliers that come under selection MORE times?
ggplot(data = con.overall[pos >= 1], aes(y = pos))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(x = 0, color = pos), width = .36)+
  scale_color_gradient(low = adjustcolor("grey80",alpha = .2), high = "blue", name = "Times under\npositive selection")+
  ylab(label = "Times Gene is positively selected")+
  labs(title = "Are there outlier genes that are under selection many more times?", 
       subtitle = "each dot is a gene, only genes under positive selection at least once are included")


# are the same genes under selection before and after?

p.pre <- ggplot(data = con.pre.post[pos.pre >= 1], aes(y = pos.pre))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(x = 0, color = pos.post), width = .36)+
  scale_color_gradient(low = adjustcolor("grey80",alpha = .8), high = "red2", name = "Times under\npositive selection\nAFTER 2012")+
  ylab(label = "Times Gene is positively selected")+
  labs(title = "Pre-2012")


p.post <- ggplot(data = con.pre.post[pos.post >= 1], aes(y = pos.post))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(x = 0, color = pos.pre), width = .36)+
  scale_color_gradient(low = adjustcolor("grey80",alpha = .8), high = "red2", name = "Times under\npositive selection\nBEFORE 2012")+
  ylab(label = "Times Gene is positively selected")+
  labs(title = "Post-2012")

p.pre + p.post + plot_annotation(title = "Are there outlier genes that are under selection many more times?", 
                                 subtitle = "each dot is a gene, only genes under positive selection at least once are included")

# most genes selected before are also selected after, but
# there are a fair number of NEW highly selected genes in the after group

# ---- excited, so just smoosh into this script the annotations, too! ----

# run 2023-09-19_convert_KO_to_pathways
genome.annotations
genome.annotations <- genome.annotations[ , .(pathway.description = paste(unique(paste0(pathway.description, " (",round(perc.complete,0),")")), collapse = "\n")), by = .(gene,ko,ko.description,threshold,score,e.value)]

# looks like some genes have more than one ko!
nrow(genome.annotations[,.N,by = .(gene)])
nrow(genome.annotations[,.N,by = .(gene,ko)]) # WTF?? oh! different significances?
genome.annotations[gene == "ME2011-09-21_3300043464_group3_bin69_scaffold_10457_c1_10",.N,by = .(gene,ko,threshold, e.value, ko.description)]

key.genes.anns <- merge(x = genome.annotations, y = key.genes.info, by = "gene", all.x = F, all.y = TRUE, allow.cartesian = TRUE)


