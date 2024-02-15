# RRR

# OK, so this is a long shot, but maybe certain KO's or at least pathways/modules are shared btwn taxa, or btwn seasons
# should at least look even though I'm not optimistic

# But... how to process. Each org is in a different table right now
# ahh, I do have the small tables of heatmap keys, could combine those into 1 not-so-big table
# need to add in the genome name to each, and then combine. 
# do this with 2023-09-25/process_KEGG_data/combine_consistent_gene_tables.R

# for the function that is called inside the data.table brackets, this explains the .SD variable: https://cran.r-project.org/web/packages/data.table/vignettes/datatable-sd-usage.html

# ---- set-up ----

library(data.table)

genes <- fread(input = "data/2023-10-16_consistent_genes_example_data/combined_file/key_genes_all_genomes-anns_separate_lines.tsv.gz")

tax <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

# ---- define functions ----

get.outliers <- function(my.vect){
  outliers <- boxplot(my.vect[ ,N.genomes], plot = F)
  outliers <- outliers$out
  if (length(outliers > 0)){
    my.vect <- my.vect[, N.genomes] >= min(outliers, na.rm = T)
  }else{
    my.vect <- rep(x = FALSE, times = nrow(my.vect))
  }
  
  # my.vect <- sum(my.vect) # having a single value output would make group by collapse rows
  return(my.vect)
}

get.key.KOs <- function(my.tab, my.tax.col, consistent.in.group, outliers.only = TRUE, subset.key.ones = TRUE){

  by.tax <- my.tab[consistent.in == consistent.in.group & is.outlier == outliers.only] 
  
  choose.cols <- c(my.tax.col, "gene", "genome", "ko", "ko.description") # data.table can take all characters, or you have to do something more complicated to make the character an internal non-character colname
  by.tax <- by.tax[ , .(ko.present = TRUE),by = choose.cols] 
  
  choose.cols <- c(my.tax.col, "genome", "ko", "ko.description")
  by.tax <- by.tax[ , .(ko.present = any(ko.present)),by = choose.cols] 
  
  choose.cols <- c(my.tax.col, "ko", "ko.description")
  by.tax <- by.tax[ , .(N.genomes = .N),by = choose.cols] 
  
  by.tax <- by.tax[ , `:=`(key.ko = get.outliers(my.vect = .SD))] 

  if(subset.key.ones){
    by.tax <- by.tax[key.ko == TRUE] 
  }
  
  index <- order(by.tax[[my.tax.col]], by.tax[["N.genomes"]], decreasing = TRUE)
  by.tax <- by.tax[index, ]

  return(by.tax)
}

get.key.modules <- function(my.tab, my.tax.col, consistent.in.group, outliers.only = TRUE, subset.key.ones = TRUE){
  
  by.tax <- my.tab[consistent.in == consistent.in.group & is.outlier == outliers.only] 
  
  choose.cols <- c(my.tax.col, "gene", "genome", "module", "module.description") # data.table can take all characters, or you have to do something more complicated to make the character an internal non-character colname
  by.tax <- by.tax[ , .(module.present = TRUE),by = choose.cols] 
  
  choose.cols <- c(my.tax.col, "genome", "module", "module.description")
  by.tax <- by.tax[ , .(module.present = any(module.present)),by = choose.cols] 
  
  choose.cols <- c(my.tax.col, "module", "module.description")
  by.tax <- by.tax[ , .(N.genomes = .N),by = choose.cols] 
  
  by.tax <- by.tax[ , `:=`(key.module = get.outliers(my.vect = .SD))] 
  
  if(subset.key.ones){
    by.tax <- by.tax[key.module == TRUE] 
  }
  
  index <- order(by.tax[[my.tax.col]], by.tax[["N.genomes"]], decreasing = TRUE)
  by.tax <- by.tax[index, ]
  
  return(by.tax)
}

get.key.pathways <- function(my.tab, my.tax.col, consistent.in.group, outliers.only = TRUE, subset.key.ones = TRUE){
  
  by.tax <- my.tab[consistent.in == consistent.in.group & is.outlier == outliers.only] 
  
  choose.cols <- c(my.tax.col, "gene", "genome", "path", "pathway.description") # data.table can take all characters, or you have to do something more complicated to make the character an internal non-character colname
  by.tax <- by.tax[ , .(pathway.present = TRUE),by = choose.cols] 
  
  choose.cols <- c(my.tax.col, "genome", "path", "pathway.description")
  by.tax <- by.tax[ , .(pathway.present = any(pathway.present)),by = choose.cols] 
  
  choose.cols <- c(my.tax.col, "path", "pathway.description")
  by.tax <- by.tax[ , .(N.genomes = .N),by = choose.cols] 
  
  by.tax <- by.tax[ , `:=`(key.path = get.outliers(my.vect = .SD))] 
  
  if(subset.key.ones){
    by.tax <- by.tax[key.path == TRUE]  
  }
  
  index <- order(by.tax[[my.tax.col]], by.tax[["N.genomes"]], decreasing = TRUE)
  by.tax <- by.tax[index, ]
  
  return(by.tax)
}

add.tax.info <- function(my.tab, tax.tab, my.tax.col){
  x <- copy(my.tab)
  
  tot <- tax.tab[ ,.(tot.genomes = .N), by = my.tax.col]
  
  x <- merge(x = x, y = tot, by = my.tax.col, all.x = TRUE, all.y = FALSE)
  
  x[ ,perc.genomes := N.genomes / tot.genomes * 100]
  
  index <- order(x[[my.tax.col]], x[["perc.genomes"]], decreasing = TRUE)
  x <- x[index, ]
  
  return(x)
}

get.season.kos <- function(my.tab, my.tax.col, tax.tab, consistent.in.group, subset.key.ones = F){
  key.KOs <- get.key.KOs(my.tab = my.tab, my.tax.col = my.tax.col, consistent.in.group = consistent.in.group, subset.key.ones = subset.key.ones)
  key.KOs <- add.tax.info(my.tab = key.KOs, tax.tab = tax.tab, my.tax.col = my.tax.col)
  key.KOs$season <- consistent.in.group
  return(key.KOs)
}

get.season.modules <- function(my.tab, my.tax.col, tax.tab, consistent.in.group, subset.key.ones = F){
  key.modules <- get.key.modules(my.tab = my.tab, my.tax.col = my.tax.col, consistent.in.group = consistent.in.group, subset.key.ones = subset.key.ones)
  key.modules <- add.tax.info(my.tab = key.modules, tax.tab = tax.tab, my.tax.col = my.tax.col)
  key.modules$season <- consistent.in.group
  return(key.modules)
}

get.season.pathways <- function(my.tab, my.tax.col, tax.tab, consistent.in.group, subset.key.ones = F){
  key.pathways <- get.key.pathways(my.tab = my.tab, my.tax.col = my.tax.col, consistent.in.group = consistent.in.group, subset.key.ones = subset.key.ones)
  key.pathways <- add.tax.info(my.tab = key.pathways, tax.tab = tax.tab, my.tax.col = my.tax.col)
  key.pathways$season <- consistent.in.group
  return(key.pathways)
}

make.excel.formatted.file <- function(my.tab, kegg.type, has.seasons = FALSE){
  for.excel <- copy(my.tab)
  for.excel[ ,stats := paste0("(",N.genomes," / ",tot.genomes," = ",round(perc.genomes, digits = 0),"%)")]
  
  if(has.seasons){
    if(kegg.type == "ko"){
      for.excel <- dcast(data = for.excel, formula = season + ko + ko.description ~ genus, value.var = "stats")
    }else if(kegg.type == "module"){
      for.excel <- dcast(data = for.excel, formula = season + module + module.description ~ genus, value.var = "stats")
    }else if(kegg.type == "pathway"){
      for.excel <- dcast(data = for.excel, formula = season + path + pathway.description ~ genus, value.var = "stats")
    }else{
      return(cat("kegg.type should be ko, module, or pathway\n"))
    }
  }else{
    if(kegg.type == "ko"){
      for.excel <- dcast(data = for.excel, formula = ko + ko.description ~ genus, value.var = "stats")
    }else if(kegg.type == "module"){
      for.excel <- dcast(data = for.excel, formula = module + module.description ~ genus, value.var = "stats")
    }else if(kegg.type == "pathway"){
      for.excel <- dcast(data = for.excel, formula = path + pathway.description ~ genus, value.var = "stats")
    }else{
      return(cat("kegg.type should be ko, module, or pathway\n"))
    }
  }
  
  for.excel <- as.matrix(for.excel)
  for.excel[is.na(for.excel)] <- ""
  x.phylum <- sub("^.*p__", "", colnames(for.excel)) 
  x.phylum <- sub(";.*$","",x.phylum)
  x.class <- sub("^.*c__", "", colnames(for.excel)) 
  x.class <- sub(";.*$","",x.class)
  x.order <- sub("^.*o__", "", colnames(for.excel)) 
  x.order <- sub(";.*$","",x.order)
  x.family <- sub("^.*f__", "", colnames(for.excel)) 
  x.family <- sub(";.*$","",x.family)
  x.genus <- sub("^.*g__", "", colnames(for.excel)) 
  x.genus <- sub(";.*$","",x.genus)
  tax.labs <- data.table("phylum" = x.phylum, "class" = x.class, "order" = x.order, "family" = x.family, "genus" = x.genus)
  tax.labs <- as.matrix(tax.labs)
  tax.labs <- t(tax.labs)
  for.excel <- rbind(tax.labs, for.excel)
  return(for.excel)
}


# ---- format data ----

# get rid of "empty" genes that occurred when empty tables for a genome were combined
genes <- genes[gene != ""]

# add tax info
tax <- as.data.table(tax)
tax <- tax[winner == TRUE]
tax <- tax[ ,c("bin.full.name","domain","phylum","class","order","family","genus","species")]
colnames(tax)

# make tax names at each level unique by concatenating with higher levels
tax[ ,`:=`(class = paste(phylum, class, sep = "; "), 
           order = paste(phylum, class, order, sep = "; "),
           family = paste(phylum, class, order, family, sep = "; "),
           genus = paste(phylum, class, order, family, genus, sep = "; "),
           species = paste(phylum, class, order, family, genus, species, sep = "; "))]

genes <- merge(x = genes, y = tax, by.x = "genome", by.y = "bin.full.name")

# get unique by KO or module or pathway (because one KO can be in multiple pathways, for example)
kos <- genes[ , .N,by = .(genome,gene,consistent.in,pres,pos,pos.perc,is.Q4,is.outlier,threshold,e.value,signif,relaxed.signif.6,
                        domain,phylum,class,order,family,genus,species,
                        ko,ko.description)]
modules <- genes[ , .N,by = .(genome,gene,consistent.in,pres,pos,pos.perc,is.Q4,is.outlier,threshold,e.value,signif,relaxed.signif.6,
                              domain,phylum,class,order,family,genus,species,
                              module,module.description)]
pathways <- genes[ , .N,by = .(genome,gene,consistent.in,pres,pos,pos.perc,is.Q4,is.outlier,threshold,e.value,signif,relaxed.signif.6,
                               domain,phylum,class,order,family,genus,species,
                               path,pathway.description)]
# Ignore unknown annotations
kos <- kos[ko != ""]
modules <- modules[module != ""]
pathways <- pathways[path != ""]

# ---- look for common ones by taxonomy ----

# # examples of how to use the functions
# key.KOs.phylum.overall <- get.key.KOs(my.tab = kos, my.tax.col = "phylum", consistent.in.group = "overall")
# 
# key.KOs.family.overall <- get.key.KOs(my.tab = kos, my.tax.col = "family", consistent.in.group = "overall")
# key.KOs.family.overall <- add.tax.info(my.tab = key.KOs.family.overall, tax.tab = tax, my.tax.col = "family")
# 
# key.KOs.genus.overall <- get.key.KOs(my.tab = kos, my.tax.col = "genus", consistent.in.group = "overall")
# key.KOs.genus.overall <- add.tax.info(my.tab = key.KOs.genus.overall, tax.tab = tax, my.tax.col = "genus")
# 
# key.modules.genus.overall <- get.key.modules(my.tab = modules, my.tax.col = "genus", consistent.in.group = "overall")
# key.modules.genus.overall <- add.tax.info(my.tab = key.modules.genus.overall, tax.tab = tax, my.tax.col = "genus")
# 
# key.paths.genus.overall <- get.key.pathways(my.tab = pathways, my.tax.col = "genus", consistent.in.group = "overall")
# key.paths.genus.overall <- add.tax.info(my.tab = key.paths.genus.overall, tax.tab = tax, my.tax.col = "genus")


# ---- look for common ones by season & genus ----
unique(kos$consistent.in)

key.KOs.genus.ice <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Ice.On")
key.KOs.genus.spring <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Spring")
key.KOs.genus.clear <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Clearwater")
key.KOs.genus.early <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Early.Summer")
key.KOs.genus.late <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Late.Summer")
key.KOs.genus.fall <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Fall")

key.KOs.genus.season.all <- rbind(key.KOs.genus.ice, key.KOs.genus.spring, key.KOs.genus.clear, key.KOs.genus.early, key.KOs.genus.late, key.KOs.genus.fall)
rm(key.KOs.genus.ice, key.KOs.genus.spring, key.KOs.genus.clear, key.KOs.genus.early, key.KOs.genus.late, key.KOs.genus.fall)

key.modules.genus.ice <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Ice.On")
key.modules.genus.spring <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Spring")
key.modules.genus.clear <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Clearwater")
key.modules.genus.early <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Early.Summer")
key.modules.genus.late <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Late.Summer")
key.modules.genus.fall <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Fall")

key.modules.genus.season.all <- rbind(key.modules.genus.ice, key.modules.genus.spring, key.modules.genus.clear, key.modules.genus.early, key.modules.genus.late, key.modules.genus.fall)
rm(key.modules.genus.ice, key.modules.genus.spring, key.modules.genus.clear, key.modules.genus.early, key.modules.genus.late, key.modules.genus.fall)

key.pathways.genus.ice <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Ice.On")
key.pathways.genus.spring <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Spring")
key.pathways.genus.clear <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Clearwater")
key.pathways.genus.early <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Early.Summer")
key.pathways.genus.late <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Late.Summer")
key.pathways.genus.fall <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Fall")

key.pathways.genus.season.all <- rbind(key.pathways.genus.ice, key.pathways.genus.spring, key.pathways.genus.clear, key.pathways.genus.early, key.pathways.genus.late, key.pathways.genus.fall)
rm(key.pathways.genus.ice, key.pathways.genus.spring, key.pathways.genus.clear, key.pathways.genus.early, key.pathways.genus.late, key.pathways.genus.fall)

# highlight only annotations in lots of genomes
key.KOs.genus.ice <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Ice.On", subset.key.ones = TRUE)
key.KOs.genus.spring <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Spring", subset.key.ones = TRUE)
key.KOs.genus.clear <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Clearwater", subset.key.ones = TRUE)
key.KOs.genus.early <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Early.Summer", subset.key.ones = TRUE)
key.KOs.genus.late <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Late.Summer", subset.key.ones = TRUE)
key.KOs.genus.fall <- get.season.kos(my.tab = kos, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Fall", subset.key.ones = TRUE)

key.KOs.genus.season.highlights <- rbind(key.KOs.genus.ice, key.KOs.genus.spring, key.KOs.genus.clear, key.KOs.genus.early, key.KOs.genus.late, key.KOs.genus.fall)
rm(key.KOs.genus.ice, key.KOs.genus.spring, key.KOs.genus.clear, key.KOs.genus.early, key.KOs.genus.late, key.KOs.genus.fall)

key.modules.genus.ice <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Ice.On", subset.key.ones = TRUE)
key.modules.genus.spring <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Spring", subset.key.ones = TRUE)
key.modules.genus.clear <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Clearwater", subset.key.ones = TRUE)
key.modules.genus.early <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Early.Summer", subset.key.ones = TRUE)
key.modules.genus.late <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Late.Summer", subset.key.ones = TRUE)
key.modules.genus.fall <- get.season.modules(my.tab = modules, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Fall", subset.key.ones = TRUE)

key.modules.genus.season.highlights <- rbind(key.modules.genus.ice, key.modules.genus.spring, key.modules.genus.clear, key.modules.genus.early, key.modules.genus.late, key.modules.genus.fall)
rm(key.modules.genus.ice, key.modules.genus.spring, key.modules.genus.clear, key.modules.genus.early, key.modules.genus.late, key.modules.genus.fall)

key.pathways.genus.ice <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Ice.On", subset.key.ones = TRUE)
key.pathways.genus.spring <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Spring", subset.key.ones = TRUE)
key.pathways.genus.clear <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Clearwater", subset.key.ones = TRUE)
key.pathways.genus.early <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Early.Summer", subset.key.ones = TRUE)
key.pathways.genus.late <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Late.Summer", subset.key.ones = TRUE)
key.pathways.genus.fall <- get.season.pathways(my.tab = pathways, my.tax.col = "genus", tax.tab = tax, consistent.in.group = "Fall", subset.key.ones = TRUE)

key.pathways.genus.season.highlights <- rbind(key.pathways.genus.ice, key.pathways.genus.spring, key.pathways.genus.clear, key.pathways.genus.early, key.pathways.genus.late, key.pathways.genus.fall)
rm(key.pathways.genus.ice, key.pathways.genus.spring, key.pathways.genus.clear, key.pathways.genus.early, key.pathways.genus.late, key.pathways.genus.fall)


# ---- export some tables for excel viewing- overall selection ----

key.KOs <- get.key.KOs(my.tab = kos, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = FALSE)
key.KOs <- add.tax.info(my.tab = key.KOs, tax.tab = tax, my.tax.col = "genus")
spreadsheet <- make.excel.formatted.file(my.tab = key.KOs, kegg.type = "ko")
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_KOs-all.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

key.KOs <- get.key.KOs(my.tab = kos, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = TRUE)
key.KOs <- add.tax.info(my.tab = key.KOs, tax.tab = tax, my.tax.col = "genus")
spreadsheet <- make.excel.formatted.file(my.tab = key.KOs, kegg.type = "ko")
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_KOs-highlights.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

key.modules <- get.key.modules(my.tab = modules, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = FALSE)
key.modules <- add.tax.info(my.tab = key.modules, tax.tab = tax, my.tax.col = "genus")
spreadsheet <- make.excel.formatted.file(my.tab = key.modules, kegg.type = "module")
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_modules-all.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

key.modules <- get.key.modules(my.tab = modules, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = TRUE)
key.modules <- add.tax.info(my.tab = key.modules, tax.tab = tax, my.tax.col = "genus")
spreadsheet <- make.excel.formatted.file(my.tab = key.modules, kegg.type = "module")
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_modules-highlights.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

key.pathways <- get.key.pathways(my.tab = pathways, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = FALSE)
key.pathways <- add.tax.info(my.tab = key.pathways, tax.tab = tax, my.tax.col = "genus")
spreadsheet <- make.excel.formatted.file(my.tab = key.pathways, kegg.type = "pathway")
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_pathways-all.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

key.pathways <- get.key.pathways(my.tab = pathways, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = TRUE)
key.pathways <- add.tax.info(my.tab = key.pathways, tax.tab = tax, my.tax.col = "genus")
spreadsheet <- make.excel.formatted.file(my.tab = key.pathways, kegg.type = "pathway")
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_pathways-highlights.tsv", sep = "\t", row.names = F, col.names = F, quote = F)


# ---- export some tables for excel viewing- seasonal selection ----

spreadsheet <- make.excel.formatted.file(my.tab = key.KOs.genus.season.all, kegg.type = "ko", has.seasons = TRUE)
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_KOs-all.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

spreadsheet <- make.excel.formatted.file(my.tab = key.KOs.genus.season.highlights, kegg.type = "ko", has.seasons = TRUE)
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_KOs-highlights.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

spreadsheet <- make.excel.formatted.file(my.tab = key.modules.genus.season.all, kegg.type = "module", has.seasons = TRUE)
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_modules-all.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

spreadsheet <- make.excel.formatted.file(my.tab = key.modules.genus.season.highlights, kegg.type = "module", has.seasons = TRUE)
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_modules-highlights.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

spreadsheet <- make.excel.formatted.file(my.tab = key.pathways.genus.season.all, kegg.type = "pathway", has.seasons = TRUE)
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_pathways-all.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

spreadsheet <- make.excel.formatted.file(my.tab = key.pathways.genus.season.highlights, kegg.type = "pathway", has.seasons = TRUE)
write.table(x = spreadsheet, file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_pathways-highlights.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

# ---- visualize overall selection ----
# overall KOs by genus (key KOs) ----

key.KOs <- get.key.KOs(my.tab = kos, my.tax.col = "genus", consistent.in.group = "overall")
key.KOs <- add.tax.info(my.tab = key.KOs, tax.tab = tax, my.tax.col = "genus")

wide.tab <- dcast(data = key.KOs, formula = ko + ko.description ~ genus, value.var = "perc.genomes")
mat <- as.matrix(wide.tab[ ,-1], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_KOs-highlights.pdf", width = 15, height = 30)

par(fig = c(0,.9,0,1), mar = c(.5,22,8,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = 1)
axis(side = 2, at = 1:ncol(image.mat), labels = colnames(image.mat), las = 2, lwd = 0, line = -.75, cex.axis = .5)


par(fig = c(.9,1,.6,.7), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
mtext(text = "Genomes with\nKO consistently\nselected (%)", adj = 0, line = 1, cex = 1, at = -.25)

dev.off()

# overall KOs by genus (all KOs) ----

key.KOs <- get.key.KOs(my.tab = kos, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = FALSE)
key.KOs <- add.tax.info(my.tab = key.KOs, tax.tab = tax, my.tax.col = "genus")

wide.tab <- dcast(data = key.KOs, formula = ko + ko.description ~ genus, value.var = "perc.genomes")
mat <- as.matrix(wide.tab[ ,-1], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_KOs-all.pdf", width = 30, height = 30)

par(fig = c(0,.94,0,1), mar = c(.5,28,8,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = 1)
axis(side = 2, at = 1:ncol(image.mat), labels = colnames(image.mat), las = 2, lwd = 0, line = -.75, cex.axis = .5)

par(fig = c(.94,.99,.6,.7), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
mtext(text = "Genomes with\nKO consistently\nselected (%)", adj = 0, line = 1, cex = 1, at = -.25)

dev.off()

# overall modules by genus (key) ----

key.modules <- get.key.modules(my.tab = modules, my.tax.col = "genus", consistent.in.group = "overall")
key.modules <- add.tax.info(my.tab = key.modules, tax.tab = tax, my.tax.col = "genus")

wide.tab <- dcast(data = key.modules, formula = module + module.description ~ genus, value.var = "perc.genomes")
mat <- as.matrix(wide.tab[ ,-1], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_modules-highlights.pdf", width = 15, height = 15)

par(fig = c(0,.9,0,1), mar = c(.5,17,8,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = 1)
axis(side = 2, at = 1:ncol(image.mat), labels = colnames(image.mat), las = 2, lwd = 0, line = -.75, cex.axis = .5)

par(fig = c(.9,1,.5,.7), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
mtext(text = "Genomes with\nmodule\nconsistently\nselected (%)", adj = 0, line = 1, cex = 1, at = -.25)

dev.off()

# overall modules by genus (all) ----

key.modules <- get.key.modules(my.tab = modules, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = FALSE)
key.modules <- add.tax.info(my.tab = key.modules, tax.tab = tax, my.tax.col = "genus")

wide.tab <- dcast(data = key.modules, formula = module + module.description ~ genus, value.var = "perc.genomes")
mat <- as.matrix(wide.tab[ ,-1], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_modules-all.pdf", width = 30, height = 30)

par(fig = c(0,.94,0,1), mar = c(.5,28,8.5,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = 1)
axis(side = 2, at = 1:ncol(image.mat), labels = colnames(image.mat), las = 2, lwd = 0, line = -.75, cex.axis = .5)

par(fig = c(.94,.99,.6,.7), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
mtext(text = "Genomes with\nmodule\nconsistently\nselected (%)", adj = 0, line = 1, cex = 1, at = -.25)

dev.off()

# overall pathways by genus (key) ----

key.pathways <- get.key.pathways(my.tab = pathways, my.tax.col = "genus", consistent.in.group = "overall")
key.pathways <- add.tax.info(my.tab = key.pathways, tax.tab = tax, my.tax.col = "genus")

wide.tab <- dcast(data = key.pathways, formula = path + pathway.description ~ genus, value.var = "perc.genomes")
mat <- as.matrix(wide.tab[ ,-1], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_pathways-highlights.pdf", width = 20, height = 10)

par(fig = c(0,.93,0,1), mar = c(.5,12,9,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = .7)
axis(side = 2, at = 1:ncol(image.mat), labels = colnames(image.mat), las = 2, lwd = 0, line = -.75, cex.axis = .7)

par(fig = c(.93,1,.5,.7), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
mtext(text = "Genomes with\npathway\nconsistently\nselected (%)", adj = 0, line = 1, cex = 1, at = -.25)

dev.off()

# overall pathways by genus (all) ----

key.pathways <- get.key.pathways(my.tab = pathways, my.tax.col = "genus", consistent.in.group = "overall", subset.key.ones = FALSE)
key.pathways <- add.tax.info(my.tab = key.pathways, tax.tab = tax, my.tax.col = "genus")

wide.tab <- dcast(data = key.pathways, formula = path + pathway.description ~ genus, value.var = "perc.genomes")
mat <- as.matrix(wide.tab[ ,-1], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/overall_genus_pathways-all.pdf", width = 30, height = 30)

par(fig = c(0,.94,0,1), mar = c(.5,12,5,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
axis(side = 2, at = 1:ncol(image.mat), labels = colnames(image.mat), las = 2, lwd = 0, line = -.75, cex.axis = .5)


par(fig = c(.94,.99,.6,.7), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
mtext(text = "Genomes with\nmodule\nconsistently\nselected (%)", adj = 0, line = 1, cex = 1, at = -.25)

dev.off()

# ---- visualize season selection ----
# season KOs by genus (key KOs) ----

col.season <- data.table("color" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), "season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))

wide.tab <- dcast(data = key.KOs.genus.season.highlights, formula = season + ko + ko.description ~ genus, value.var = "perc.genomes")
wide.tab <- merge(x = wide.tab, y = col.season, by = "season")
wide.tab[ ,season := factor(x = season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), ordered = T)]
wide.tab <- wide.tab[order(season)]
mat.columns <- 3:(ncol(wide.tab) - 1)
mat <- as.matrix(wide.tab[ ,..mat.columns], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))
kegg.labs <- colnames(image.mat)
all.equal(kegg.labs, wide.tab$ko.description[nrow(wide.tab):1])
kegg.colors <- wide.tab$color[nrow(wide.tab):1]

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_KOs-highlights.pdf", width = 11, height = 15)

par(fig = c(0,.88,0,1), mar = c(.5,28,5,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
axis(side = 2, at = 1:ncol(image.mat), labels = kegg.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
index <- which(kegg.colors == "snow3")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "snow3")
index <- which(kegg.colors == "tan4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "tan4")
index <- which(kegg.colors == "cornflowerblue")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "cornflowerblue")
index <- which(kegg.colors == "chartreuse4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "chartreuse4")
index <- which(kegg.colors == "purple")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "purple")
index <- which(kegg.colors == "hotpink2")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "hotpink2")

par(fig = c(.88,1,.6,.7), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
mtext(text = "Genomes with\nKO consistently\nselected (%)", adj = 0, line = 1, cex = .9, at = -.5)

dev.off()

# season KOs by genus (all KOs) ----

col.season <- data.table("color" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), "season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))

wide.tab <- dcast(data = key.KOs.genus.season.all, formula = season + ko + ko.description ~ genus, value.var = "perc.genomes")
wide.tab <- merge(x = wide.tab, y = col.season, by = "season")
wide.tab[ ,season := factor(x = season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), ordered = T)]
wide.tab <- wide.tab[order(season)]
mat.columns <- 3:(ncol(wide.tab) - 1)
mat <- as.matrix(wide.tab[ ,..mat.columns], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))
kegg.labs <- colnames(image.mat)
all.equal(kegg.labs, wide.tab$ko.description[nrow(wide.tab):1])
kegg.colors <- wide.tab$color[nrow(wide.tab):1]

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_KOs-all.pdf", width = 30, height = 240)

par(fig = c(0,.92,0,1), mar = c(.5,25,5,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
axis(side = 2, at = 1:ncol(image.mat), labels = kegg.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
index <- which(kegg.colors == "snow3")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "snow3")
index <- which(kegg.colors == "tan4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "tan4")
index <- which(kegg.colors == "cornflowerblue")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "cornflowerblue")
index <- which(kegg.colors == "chartreuse4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "chartreuse4")
index <- which(kegg.colors == "purple")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "purple")
index <- which(kegg.colors == "hotpink2")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "hotpink2")

par(fig = c(.92,.99,.8,.82), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = 1)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F, cex.axis = .7)
mtext(text = "Genomes with\nKO consistently\nselected (%)", adj = 0, line = 1, cex = .9, at = 0)

dev.off()

# season modules by genus (key modules) ----

col.season <- data.table("color" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), "season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))

wide.tab <- dcast(data = key.modules.genus.season.highlights, formula = season + module + module.description ~ genus, value.var = "perc.genomes")
wide.tab <- merge(x = wide.tab, y = col.season, by = "season")
wide.tab[ ,season := factor(x = season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), ordered = T)]
wide.tab <- wide.tab[order(season)]
mat.columns <- 3:(ncol(wide.tab) - 1)
mat <- as.matrix(wide.tab[ ,..mat.columns], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))
kegg.labs <- colnames(image.mat)
all.equal(kegg.labs, wide.tab$ko.description[nrow(wide.tab):1])
kegg.colors <- wide.tab$color[nrow(wide.tab):1]

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_modules-highlights.pdf", width = 9, height = 18)

par(fig = c(0,.89,0,1), mar = c(.5,17,5,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
axis(side = 2, at = 1:ncol(image.mat), labels = kegg.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
index <- which(kegg.colors == "snow3")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "snow3")
index <- which(kegg.colors == "tan4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "tan4")
index <- which(kegg.colors == "cornflowerblue")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "cornflowerblue")
index <- which(kegg.colors == "chartreuse4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "chartreuse4")
index <- which(kegg.colors == "purple")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "purple")
index <- which(kegg.colors == "hotpink2")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "hotpink2")

par(fig = c(.87,1,.6,.7), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = .6)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
mtext(text = "Genomes with\nmodule\nconsistently\nselected (%)", adj = 0, line = 1, cex = .6, at = 0)

dev.off()

# season modules by genus (all modules) ----

col.season <- data.table("color" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), "season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))

wide.tab <- dcast(data = key.modules.genus.season.all, formula = season + module + module.description ~ genus, value.var = "perc.genomes")
wide.tab <- merge(x = wide.tab, y = col.season, by = "season")
wide.tab[ ,season := factor(x = season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), ordered = T)]
wide.tab <- wide.tab[order(season)]
mat.columns <- 3:(ncol(wide.tab) - 1)
mat <- as.matrix(wide.tab[ ,..mat.columns], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))
kegg.labs <- colnames(image.mat)
all.equal(kegg.labs, wide.tab$ko.description[nrow(wide.tab):1])
kegg.colors <- wide.tab$color[nrow(wide.tab):1]

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_modules-all.pdf", width = 20, height = 90)

par(fig = c(0,.93,0,1), mar = c(.5,17,5,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
axis(side = 2, at = 1:ncol(image.mat), labels = kegg.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
index <- which(kegg.colors == "snow3")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "snow3")
index <- which(kegg.colors == "tan4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "tan4")
index <- which(kegg.colors == "cornflowerblue")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "cornflowerblue")
index <- which(kegg.colors == "chartreuse4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "chartreuse4")
index <- which(kegg.colors == "purple")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "purple")
index <- which(kegg.colors == "hotpink2")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "hotpink2")

par(fig = c(.93,1,.8,.83), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = .6)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
mtext(text = "Genomes with\nmodule\nconsistently\nselected (%)", adj = 0, line = 1, cex = .6, at = 0)

dev.off()

# season pathways by genus (key pathways) ----

col.season <- data.table("color" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), "season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))

wide.tab <- dcast(data = key.pathways.genus.season.highlights, formula = season + path + pathway.description ~ genus, value.var = "perc.genomes")
wide.tab <- merge(x = wide.tab, y = col.season, by = "season")
wide.tab[ ,season := factor(x = season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), ordered = T)]
wide.tab <- wide.tab[order(season)]
mat.columns <- 3:(ncol(wide.tab) - 1)
mat <- as.matrix(wide.tab[ ,..mat.columns], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))
kegg.labs <- colnames(image.mat)
all.equal(kegg.labs, wide.tab$ko.description[nrow(wide.tab):1])
kegg.colors <- wide.tab$color[nrow(wide.tab):1]

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_pathways-highlights.pdf", width = 15, height = 60)

par(fig = c(0,.91,0,1), mar = c(.5,12,5,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
axis(side = 2, at = 1:ncol(image.mat), labels = kegg.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
index <- which(kegg.colors == "snow3")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "snow3")
index <- which(kegg.colors == "tan4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "tan4")
index <- which(kegg.colors == "cornflowerblue")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "cornflowerblue")
index <- which(kegg.colors == "chartreuse4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "chartreuse4")
index <- which(kegg.colors == "purple")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "purple")
index <- which(kegg.colors == "hotpink2")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "hotpink2")

par(fig = c(.91,1,.8,.85), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = .6)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
mtext(text = "Genomes with\nmodule\nconsistently\nselected (%)", adj = 0, line = 1, cex = .6, at = 0)

dev.off()

# season pathways by genus (all pathways) ----

col.season <- data.table("color" = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), "season" = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))

wide.tab <- dcast(data = key.pathways.genus.season.all, formula = season + path + pathway.description ~ genus, value.var = "perc.genomes")
wide.tab <- merge(x = wide.tab, y = col.season, by = "season")
wide.tab[ ,season := factor(x = season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"), ordered = T)]
wide.tab <- wide.tab[order(season)]
mat.columns <- 3:(ncol(wide.tab) - 1)
mat <- as.matrix(wide.tab[ ,..mat.columns], rownames = TRUE)

n.colors <- 100
color.fun <- colorRampPalette(colors = c("grey", "blue"))
col.key <- data.table("color" = color.fun(n = n.colors), "value" = seq(from = min(mat, na.rm = T), to = max(mat, na.rm = T), along.with = 1:n.colors))

image.mat <- t(mat)
image.mat <- image.mat[ ,ncol(image.mat):1, drop = F] # this way it matches the table key 
tax.labs <- sub(pattern = "^.*__", replacement = "", x = rownames(image.mat))
kegg.labs <- colnames(image.mat)
all.equal(kegg.labs, wide.tab$ko.description[nrow(wide.tab):1])
kegg.colors <- wide.tab$color[nrow(wide.tab):1]

pdf(file = "figures/2023-10-16_cross-genome_consistently_selected_genes_summary_attempt/season_genus_pathways-all.pdf", width = 30, height = 120)

par(fig = c(0,.95,0,1), mar = c(.5,12,5,0))
image(z = image.mat, x = 1:nrow(image.mat), y = 1:ncol(image.mat), axes = F, ann = F, col = col.key$color)
box(which = "plot", lwd = 1)
axis(side = 3, at = 1:nrow(image.mat), labels = tax.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
axis(side = 2, at = 1:ncol(image.mat), labels = kegg.labs, las = 2, lwd = 0, line = -.75, cex.axis = .5)
index <- which(kegg.colors == "snow3")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "snow3")
index <- which(kegg.colors == "tan4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "tan4")
index <- which(kegg.colors == "cornflowerblue")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "cornflowerblue")
index <- which(kegg.colors == "chartreuse4")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "chartreuse4")
index <- which(kegg.colors == "purple")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "purple")
index <- which(kegg.colors == "hotpink2")
rect(ybottom = min((1:ncol(image.mat))[index]) - .45, ytop = max((1:ncol(image.mat))[index]) + .45, xleft = 1 - .5, xright = nrow(image.mat) + .5, lwd = 4, border = "hotpink2")

par(fig = c(.95,1,.9,.93), mar = c(1,2,1,3), new = T)
plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, las = 2, line = 0, cex.axis = .6)
axis(side = 4, at = pretty(x = col.key$value), xpd = T, lwd = 0, lwd.ticks = 1, tck = -.25, line = 0, labels = F)
mtext(text = "Genomes with\nmodule\nconsistently\nselected (%)", adj = 0, line = 1, cex = .6, at = 0)

dev.off()
