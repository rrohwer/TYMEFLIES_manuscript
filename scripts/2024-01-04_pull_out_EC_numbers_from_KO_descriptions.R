# RRR

library(data.table)

key <- fread("figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/consistently_selected_genes_KEY-ME2011-09-21_3300043464_group3_bin69.csv")

key$ko.description # this has the EC numbers

key[ ,EC.number := sub(pattern = "^.*\\[EC", replacement = "EC", x = ko.description)]
key[ ,EC.number := sub(pattern = "\\]$", replacement = "", x = EC.number)]
index <- grep(pattern = "^EC", x = key$EC.number, value = F, invert = T)
key[index, EC.number := ""]
key[ ,EC.number := sub(pattern = "^.*\\EC:", replacement = "", x = EC.number)]


fwrite(x = key, file = "data/2023-12-18_manual_acI-B_metabolism_checking/acI-B_gene_annotation_files/acI-B_heatmap_key_with_EC_numbers.csv")



genes <- fread("data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.ko-pathways_anns-separate.tsv.gz")

genes[ ,EC.number := sub(pattern = "^.*\\[EC", replacement = "EC", x = ko.description)]
genes[ ,EC.number := sub(pattern = "\\]$", replacement = "", x = EC.number)]
index <- grep(pattern = "^EC", x = genes$EC.number, value = F, invert = T)
genes[index, EC.number := ""]
genes[ ,EC.number := sub(pattern = "^.*\\EC:", replacement = "", x = EC.number)]

fwrite(x = genes, file = "data/2023-12-18_manual_acI-B_metabolism_checking/acI-B_gene_annotation_files/acI-B_all_annotations_with_EC_numbers.csv")
