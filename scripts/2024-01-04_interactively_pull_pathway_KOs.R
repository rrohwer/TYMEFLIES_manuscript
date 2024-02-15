# RRR
# start here: https://www.genome.jp/brite/query=03420&htext=br08901.keg&option=-a&node_proc=br08901_org&proc_enabled=ko&panel=collapse
# enter pathway, without the "map" into the "ID search" box
# click the red link that appears on the right
# then enter the KO number, this time with the "K" in the ID search box to highlight the selected step
# then paste in all the space-separated KOs in the generated temp file to see if steps surrounding the selected step are in the genome.
# In this way I tried to figure out which pathway/big-picture process was most likely for each selected gene. 
# Because the selected gene's enzyme should lie in an existing genomic pathway. 


library(data.table)

key <- fread(file = "data/2023-12-18_manual_acI-B_metabolism_checking/acI-B_gene_annotation_files/acI-B_heatmap_key_with_EC_numbers.csv")
genes <- fread(file = "data/2023-12-18_manual_acI-B_metabolism_checking/acI-B_gene_annotation_files/acI-B_all_annotations_with_EC_numbers.csv")

temp <- genes[path == "map00020", ko]
temp <- paste(temp, collapse = " ")

write.table(x = temp, file = "data/2023-12-18_manual_acI-B_metabolism_checking/temp.txt", sep = "\n", quote = F, row.names = F, col.names = F)

