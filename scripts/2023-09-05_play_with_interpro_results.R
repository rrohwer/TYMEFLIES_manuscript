# RRR

# look at the interpro results of the 3 acI cherries
# what are the genes under selection? maybe make a table of just them.

# ---- set-up ----

library(data.table)

# B
interpro <- fread(file = "data/2023-09-05_example_interpro_data/ME2011-09-21_3300043464_group3_bin69.tsv", quote = "", sep = "\t", fill = T)
genes <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz")

# C
interpro <- fread(file = "data/2023-09-05_example_interpro_data/ME2016-07-20_3300033996_group7_bin32.tsv", quote = "", sep = "\t", fill = T)
genes <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2016-07-20_3300033996_group7_bin32_gene_info_combined_MK.tsv.gz")

# A
interpro <- fread(file = "data/2023-09-05_example_interpro_data/ME2011-09-04_3300044729_group3_bin142.tsv", quote = "", sep = "\t", fill = T)
genes <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-04_3300044729_group3_bin142_gene_info_combined_MK.tsv.gz")

# ---- get interpro data ----

colnames(interpro) <- c("Protein.Accession","MD5","Sequence.Length","Analysis","Signature.Accession","Signature.Description","Start.Location",
                        "Stop.Location","e.value","Is.Match","Run.Date","Interpro.Accession","Interpro.Description","GO.Annotations","Pathways.Annotations")

unique(interpro$Is.Match)
interpro <- interpro[ ,-c("MD5","Run.Date")]
interpro <- interpro[ ,c("Protein.Accession","Sequence.Length","Start.Location","Stop.Location","e.value","Is.Match","Analysis","Signature.Accession","Signature.Description","Interpro.Accession","Interpro.Description","GO.Annotations","Pathways.Annotations")]

# ---- get genes under selection data ----

colnames(genes)

genes <- genes[ ,c("gene","sample","date","year","yday","season","invasion",
                   "mcdonald.kreitman","MK.p.val","gene_length", "start","end","coverage","breadth","SNV_count","SNS_count","divergent_site_count",
                   "Description","COG_Description","Broad_category")]

# ---- play with it? ----

big <- merge(x = genes, y = interpro, by.x = "gene", by.y = "Protein.Accession", all = T, allow.cartesian = TRUE) # interpro has multiple annotations per gene so this blows up the table

pos <- big[mcdonald.kreitman < 1 & MK.p.val <= .05 & !is.na(mcdonald.kreitman) & is.finite(mcdonald.kreitman)]
neg <- big[mcdonald.kreitman > 1 & MK.p.val <= .05 & !is.na(mcdonald.kreitman) & is.finite(mcdonald.kreitman)]

# make a table to look at...so in order by gene but annotations only listed once

pos <- pos[ ,.(.N, date = paste(date, collapse = " "), season = paste(season, collapse = " ")),by = .(gene,Sequence.Length,Start.Location,Stop.Location,Description,COG_Description,Broad_category,Analysis,Signature.Accession,Signature.Description,Interpro.Accession,Interpro.Description,GO.Annotations,Pathways.Annotations)]
neg <- neg[ ,.(.N, date = paste(date, collapse = " "), season = paste(season, collapse = " ")),by = .(gene,Sequence.Length,Start.Location,Stop.Location,Description,COG_Description,Broad_category,Analysis,Signature.Accession,Signature.Description,Interpro.Accession,Interpro.Description,GO.Annotations,Pathways.Annotations)]

pos <- pos[ ,c("gene","Sequence.Length","Start.Location","Stop.Location","N","date","season","Description","COG_Description","Broad_category","Analysis","Signature.Accession","Signature.Description","Interpro.Accession","Interpro.Description","GO.Annotations","Pathways.Annotations")]
neg <- neg[ ,c("gene","Sequence.Length","Start.Location","Stop.Location","N","date","season","Description","COG_Description","Broad_category","Analysis","Signature.Accession","Signature.Description","Interpro.Accession","Interpro.Description","GO.Annotations","Pathways.Annotations")]

# ---- look at these in excel??? ----

# B
fwrite(x = pos, file = "data/2023-09-05_example_interpro_data/under_selection_tables/acI-B_positive.tsv", sep = "\t", quote = F, row.names = F)
fwrite(x = neg, file = "data/2023-09-05_example_interpro_data/under_selection_tables/acI-B_purifying.tsv", sep = "\t", quote = F, row.names = F)

# C
fwrite(x = pos, file = "data/2023-09-05_example_interpro_data/under_selection_tables/acI-C_positive.tsv", sep = "\t", quote = F, row.names = F)
fwrite(x = neg, file = "data/2023-09-05_example_interpro_data/under_selection_tables/acI-C_purifying.tsv", sep = "\t", quote = F, row.names = F)

# A
fwrite(x = pos, file = "data/2023-09-05_example_interpro_data/under_selection_tables/acI-A_positive.tsv", sep = "\t", quote = F, row.names = F)
fwrite(x = neg, file = "data/2023-09-05_example_interpro_data/under_selection_tables/acI-A_purifying.tsv", sep = "\t", quote = F, row.names = F)
