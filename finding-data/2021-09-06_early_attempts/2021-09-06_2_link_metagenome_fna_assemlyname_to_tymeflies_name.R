# RRR

stupid <- read.csv(file = "finding-data/2021-09-06_1_metagenome_fna_filenames_from_neha.txt", header = F)
str(stupid)
stupid[1:5, ]

key <- readRDS(file = "data/2021-05-24_metaG_stats/all_img_assemblies.rds")
# looks like these filenames are taxon IDs

colnames(stupid) <- "Filename.From.Neha"
stupid[1:5,]
sub(pattern = "\\.assembled\\.fna", replacement = "", x = stupid$Filename.From.Neha)
stupid <- cbind(stupid, "Taxon.ID" = sub(pattern = "\\.assembled\\.fna", replacement = "", x = stupid$Filename.From.Neha))

str(stupid)
key <- key[ ,c(1,8)]
str(key)
key <- key[!is.na(key$IMG.Taxon.ID), ]
key <- unique(key)

sl.less.stupid <- merge(x = stupid, y = key, by.x = "Taxon.ID", by.y = "IMG.Taxon.ID", all.x = TRUE, all.y = FALSE)

str(sl.less.stupid)


# # solved this by making key unique
# wtf <- sl.less.stupid[duplicated(sl.less.stupid$Taxon.ID), ]
# 
# index <- duplicated(c(wtf$Taxon.ID, sl.less.stupid$Taxon.ID))
# index <- index[-(1:nrow(wtf))]
# 
# wtaf1 <- sl.less.stupid[index, ]
# 
# 
# wtf <- sl.less.stupid[duplicated(sl.less.stupid$Filename.From.Neha), ]
# 
# index <- duplicated(c(wtf$Filename.From.Neha, sl.less.stupid$Filename.From.Neha))
# index <- index[-(1:nrow(wtf))]
# 
# wtaf2 <- sl.less.stupid[index, ]
# 
# 
# wtf <- sl.less.stupid[duplicated(sl.less.stupid$TYMEFLIES.Name), ]
# 
# index <- duplicated(c(wtf$TYMEFLIES.Name, sl.less.stupid$TYMEFLIES.Name))
# index <- index[-(1:nrow(wtf))]
# 
# wtaf3 <- sl.less.stupid[index, ]
# 
# sum(duplicated(stupid$Filename.From.Neha))




sl.less.stupid <- sl.less.stupid[ ,c(2,1,3)]

sl.less.stupid <- cbind(sl.less.stupid, "My.Filename" = paste0(sl.less.stupid$TYMEFLIES.Name, ".fasta"))

# write.csv(x = sl.less.stupid, file = "finding-data/2021-09-06_3_filename_key_nehas_assemblies.csv", quote = F, row.names = F)
