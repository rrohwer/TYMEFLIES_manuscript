data.folder <- "data/2022-09-22_gtdb-tk_output/gtdbtk_bac120_summary_files/"

combo <- NULL
for(f in list.files(data.folder)){
  my.tab <- read.delim(file = file.path(data.folder,f))
  combo <- rbind(combo, my.tab)
}

colnames(combo)

bin.ID <- combo[ ,1]
tax <- combo[ ,2]

taxonomy <- matrix(data = NA, nrow = nrow(combo), ncol = 8)
colnames(taxonomy) <- c("binID","domain","phylum","class","order","family","genus","species")
taxonomy[ ,1] <- bin.ID

tax <- strsplit(x = tax, split = ";")
for (r in 1:length(tax)){
  taxonomy[r,2] <- tax[[r]][1]
  taxonomy[r,3] <- tax[[r]][2]
  taxonomy[r,4] <- tax[[r]][3]
  taxonomy[r,5] <- tax[[r]][4]
  taxonomy[r,6] <- tax[[r]][5]
  taxonomy[r,7] <- tax[[r]][6]
  taxonomy[r,8] <- tax[[r]][7]
}

head(taxonomy)

saveRDS(object = taxonomy, file = "data/2022-09-22_gtdb-tk_output/bin_taxonomy.rds")
