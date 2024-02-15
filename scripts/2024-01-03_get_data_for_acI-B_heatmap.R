# RRR
# this may be a bit of a placeholder...
# made a combined all-genome table

# pull from server: consistent_genes_analysis/key_genes-by_sample/ME2011-09-21_3300043464_group3_bin69.consistently_selected_genes-by_sample.tsv.gz

library(data.table)
library(lubridate)

genes <- fread("data/2023-12-29_genome_examples_data_for_paper_figures/ME2011-09-21_3300043464_group3_bin69.consistently_selected_genes-by_sample.tsv.gz", colClasses = c("date" = "character"))

# ---- pull out data ----

genes[ ,date := parse_date_time(date, "ymd")]

genes <- genes[consistent.in == "pre-2012" | consistent.in == "post-2012"]

genes[!is.finite(mcdonald.kreitman), mcdonald.kreitman := NA] 
max(genes$mcdonald.kreitman, na.rm = T)
min(genes$mcdonald.kreitman, na.rm = T)
# so both positive a negative selection in the table, but only looked at positive to determine over-representation

# only call it pos selection (in figure) on a given day if p-value is <= .25, otherwise p-value is 1
genes$positive.selection.pvalue <- as.numeric(1)
genes[mcdonald.kreitman < 1 & MK.p.val <= .25, positive.selection.pvalue := MK.p.val]

# subset to outliers only (not just Q4)
# genes <- genes[is.outlier == TRUE]

# make a matrix
genes <- genes[order(date)]

genes <- dcast(data = genes, formula = gene + consistent.in + ko + ko.description + module + module.description + path + pathway.description ~ sample, value.var = "positive.selection.pvalue")

# put in order
genes[ ,consistent.in := factor(consistent.in, levels = c("pre-2012","post-2012"), ordered = T)]
genes <- genes[order(consistent.in,path,module,ko)]

# subset to only the key
genes.key <- genes[ ,.(consistent.in = paste(consistent.in, collapse = "; ")), by = .(gene,ko,ko.description,module,module.description,path,pathway.description)]

# if genes were selected before 2012 too, then call them overall. don't use "overall" category b/c 2012 is so strong, all those genes end up overall selected too.
genes.key <- genes.key[consistent.in != "post-2012" & consistent.in != "pre-2012", consistent.in := "overall"]
genes.key <- genes.key[ ,consistent.in := factor(consistent.in, levels = c("overall","pre-2012","post-2012"), ordered = T)]
genes.key <- genes.key[order(consistent.in, path, module, ko)]

# remove duplicate genes, as these should not be included as post-2012
genes <- genes[ ,-c("ko","ko.description","module","module.description","path","pathway.description")]
genes$consistent.in # in order with pre-2012 first
genes <- genes[duplicated(gene) == FALSE]

# put table and key in the key's order
genes.key$row.order <- 1:nrow(genes.key)
genes <- merge(x = genes.key[ ,.(gene,row.order)], y = genes, by = "gene")
genes <- genes[order(row.order)]

# make matrix for plotting
genes <- as.matrix(genes[ ,-c("row.order","consistent.in")], rownames = T)
genes[,1:5]

# go back from p-val = 1 to high p-vals are NA
genes[genes == 1] <- NA

# check order one last time
all.equal(rownames(genes), genes.key$gene) # good

# check order of sample dates
colnames(genes)

# ---- save data for plotting ----

saveRDS(object = genes, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/consistently_selected_genes-ME2011-09-21_3300043464_group3_bin69.rds")
fwrite(x = genes.key, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/consistently_selected_genes_KEY-ME2011-09-21_3300043464_group3_bin69.csv")

# add manual annotations to the genes.key table before plotting!
