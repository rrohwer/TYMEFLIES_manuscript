# RRR
# Get which pathways and modules are in the genome.
# For every KO, add all pathway and module options along with their size & completeness

# Decisions:
# 1. What is a significant KO annotation?
#       If use the KEGG star (default kofamscan significance) there are basically no complete pathways in the acI
#       Katy and Pedro use e-value as a secondary cutoff (1e-5)
#       Anvi'o relaxes the threshold for bit scores to include. https://anvio.org/help/main/programs/anvi-run-kegg-kofams/
#       Looking at the defaults, the threshold is a manual decision for each hmmer family about how good it should match (how good a bit score)
#       And the default e-value is high, just < .01
#       So defaults are relying more on threshold/bit score than e-value to say it's significant. 
#       So I will do that too. keep e-value cutoff at < .01 and give second-tier classifications by relaxing threshold.
#       I chose .6 for my threshold relaxation- that's actually more stringent than anvi'o's .5 choice.
# 2. What is a complete pathway?
#       Partly depends on how many enzymes are in the pathway!
#       I'm going to set this aside for now. Will revisit when I choose my pathway key to use.
# 3. Which pathways reference to use? 
#       For now, just using the default kegg-defined pathways. 
#       In future, could sub in other people's specific pathways instead, like the keggdecoder ones.

# example terminal call:
# Rscript add_pathway_info_to_KO_annotations.R KEGG_lists/ko_pathway_links.txt KEGG_lists/ko_module_links.txt KEGG_lists/module_pathway_links.txt KEGG_lists/pathway_list.txt KEGG_lists/ko_list.txt KEGG_lists/module_list.txt ../yggshare/current_projects/TYMEFLIES/tymeflies/kegg_annotations/kofamscan_kegg_annotations/ME2011-09-21_3300043464_group3_bin69.kofamscan.tsv.gz kegg_annotations_with_pathways-one_line_per_gene/ME2011-09-21_3300043464_group3_bin69.ko-pathways_genes-unique.tsv.gz kegg_annotations_with_pathways/ME2011-09-21_3300043464_group3_bin69.ko-pathways_anns-separate.tsv.gz 1

# ---- set-up ----

library(data.table)

userprefs <- commandArgs(trailingOnly = TRUE)
pathways.key <- userprefs[1]
module.ko.key <- userprefs[2]
module.path.key <- userprefs[3]
pathways.names <- userprefs[4]
ko.names <- userprefs[5]
module.names <- userprefs[6]
genome.annotations <- userprefs[7]
output.annotations.file.one.gene.per.line <- userprefs[8]
output.annotations.file.one.ann.per.line <- userprefs[9]
num.threads <- as.numeric(userprefs[10])
# defined in script: relax KO bit score to 60% of threshold, only keep top annotation per gene

cat("\nprocessing ",genome.annotations,"\n")

# # local test-paths
# pathways.key <- "data/2023-09-13_KEGG_annotations/KEGG_lists/ko_pathway_links.txt"
# module.ko.key <- "data/2023-09-13_KEGG_annotations/KEGG_lists/ko_module_links.txt"
# module.path.key <- "data/2023-09-13_KEGG_annotations/KEGG_lists/module_pathway_links.txt"
# pathways.names <- "data/2023-09-13_KEGG_annotations/KEGG_lists/pathway_list.txt"
# ko.names <- "data/2023-09-13_KEGG_annotations/KEGG_lists/ko_list.txt"
# module.names <- "data/2023-09-13_KEGG_annotations/KEGG_lists/module_list.txt"
# genome.annotations <- "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.kofamscan.tsv"
# output.annotations.file.one.gene.per.line <- "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.ko-pathways_genes-unique.tsv.gz"
# output.annotations.file.one.ann.per.line <- "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.ko-pathways_anns-separate.tsv.gz"
# num.threads <- 1
# 
# genome.annotations <- "data/2023-09-13_KEGG_annotations/example_files/ME2015-09-13_3300035666_group6_bin236.kofamscan.tsv.gz"

# pathways.key <-  "KEGG_lists/ko_pathway_links.txt"
# module.ko.key <- "KEGG_lists/ko_module_links.txt"
# module.path.key <- "KEGG_lists/module_pathway_links.txt"
# pathways.names <- "KEGG_lists/pathway_list.txt"
# ko.names <- "KEGG_lists/ko_list.txt"
# module.names <- "KEGG_lists/module_list.txt"
# genome.annotations <- "../yggshare/current_projects/TYMEFLIES/tymeflies/kegg_annotations/kofamscan_kegg_annotations/ME2012-05-05_3300042888_group4_bin209.kofamscan.tsv.gz"
# output.annotations.file.one.gene.per.line <- "kegg_annotations_with_pathways-one_line_per_gene/ME2012-05-05_3300042888_group4_bin209.ko-pathways_genes-unique.tsv.gz"
# output.annotations.file.one.ann.per.line <- "kegg_annotations_with_pathways/ME2012-05-05_3300042888_group4_bin209.ko-pathways_anns-separate.tsv.gz"
# num.threads <- 1

# import files
pathways.key <- fread(file = pathways.key, header = F, sep = "\t", nThread = num.threads)
module.ko.key <- fread(file = module.ko.key, header = F, sep = "\t", nThread = num.threads)
module.path.key <- fread(file = module.path.key, header = F, sep = "\t", nThread = num.threads)
pathways.names <- fread(file = pathways.names, header = F, sep = "\t", nThread = num.threads)
ko.names <- fread(file = ko.names, header = F, sep = "\t", nThread = num.threads)
module.names <- fread(file = module.names, header = F, sep = "\t", nThread = num.threads)
genome.annotations <- fread(file = genome.annotations, header = T, sep = "\t", nThread = num.threads)


# ---- format files ----

colnames(pathways.key) <- c("ko","path")
pathways.key[ ,`:=`(ko = sub(pattern = "ko:", replacement = "", x = ko),
                    path = sub(pattern = "path:", replacement = "", x = path))]
pathways.key <- pathways.key[grep(pattern = "map", x = path)]

colnames(module.ko.key) <- c("ko","module")
module.ko.key[ ,`:=`(ko = sub(pattern = "ko:", replacement = "", x = ko),
                     module = sub(pattern = "md:", replacement = "", x = module))]

colnames(module.path.key) <- c("path","module")
module.path.key[ ,`:=`(path = sub(pattern = "path:", replacement = "", x = path),
                       module = sub(pattern = "md:", replacement = "", x = module))]

colnames(pathways.names) <- c("path","pathway.description")

colnames(ko.names) <- c("ko","ko.description")

colnames(module.names) <- c("module","module.description")

genome.annotations <- genome.annotations[-1]
colnames(genome.annotations) <- c("signif","gene","ko","threshold","bit.score","e.value","ko.description")
genome.annotations[ ,`:=`(threshold = as.numeric(threshold), bit.score = as.numeric(bit.score), e.value = as.numeric(e.value))]

# make a relaxed significance cutoff- 60% of the default threshold ----
genome.annotations$relaxed.signif.6 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .6), relaxed.signif.6 := TRUE]

# choose only best KO per gene, no multiple annotations remaining ----
genome.annotations <- genome.annotations[relaxed.signif.6 == TRUE]

genome.annotations[ ,perc.threshold := (bit.score - threshold) / threshold] # choose the best based on how close the score was to the threshold, normalized by the threshold
choose.best <- genome.annotations[ , .(perc.threshold = max(perc.threshold)),by = .(gene)]
if (any(duplicated(choose.best$gene))){cat("\n\nOH NO DUPLICATE GENES IN THE TABLE\n\n")}

genome.annotations <- merge(x = choose.best, y = genome.annotations, by = c("gene","perc.threshold"), all.x = T, all.y = F)
if (any(duplicated(genome.annotations$gene))){ # this happens if there are TIES in the perc.threshold for 2 annotations (happened twice in all genomes)
  index <- which(genome.annotations$gene == genome.annotations$gene[duplicated(genome.annotations$gene)])
  keep.index <- which(genome.annotations$e.value[index] == min(genome.annotations$e.value[index])) # choose the lower evalue one in case of bit score-based tie
  remove.index <- index[-keep.index]
  genome.annotations <- genome.annotations[-remove.index,]
} 
if (any(duplicated(genome.annotations$gene))){cat("\n\nOH NO DUPLICATE GENES IN THE TABLE\n\n")}


# ---- which pathways are in the genome ----

# for each pathway, list its KOs 

pathways.list <- list(NULL)
for (p in unique(pathways.key$path)){
  k <- pathways.key[path == p]
  k <- list(k$ko)
  names(k) <- p
  pathways.list <- c(pathways.list, k)
}

# for each pathway, how many of its KOs are present?

genome.pathways <- pathways.key[ , .N, by = .(path)]
genome.pathways$genome.N <- 0

for (p in genome.pathways$path){
  path.kos <- pathways.list[[p]]
  x <- intersect(x = path.kos, genome.annotations$ko)
  genome.pathways[path == p, genome.N := length(x)]
}

genome.pathways[ ,perc.complete := genome.N / N * 100]

# ---- which modules are in the genome ----

# for each module, list its KOs 

modules.list <- list(NULL)
for (m in unique(module.ko.key$module)){
  k <- module.ko.key[module == m]
  k <- list(k$ko)
  names(k) <- m
  modules.list <- c(modules.list, k)
}

# for each module, how many of its KOs are present?

genome.modules <- module.ko.key[ , .N, by = .(module)]
genome.modules$genome.N <- 0

for (m in genome.modules$module){
  module.kos <- modules.list[[m]]
  x <- intersect(x = module.kos, genome.annotations$ko)
  genome.modules[module == m, genome.N := length(x)]
}

genome.modules[ ,perc.complete := genome.N / N * 100]


# put names to the pathways and modules and sort by how big they are ----

genome.pathways <- merge(x = pathways.names, y = genome.pathways, by = "path", all.x = F, all.y = T)

genome.pathways <- genome.pathways[order(N, decreasing = T)]

genome.modules <- merge(x = module.names, y = genome.modules, by = "module", all.x = F, all.y = T)

genome.modules <- genome.modules[order(N, decreasing = T)]


# ---- which pathways do the genome's KOs belong to ----

# only include genes with known annotation, but include even if no assigned pathway
genome.annotations <- merge(x = genome.annotations, y = pathways.key, by = "ko", all.x = T, all.y = F, allow.cartesian = TRUE)
genome.annotations <- merge(x = genome.annotations, y = genome.pathways, by = "path", all.x = T, all.y = F)
colnames(genome.annotations)[colnames(genome.annotations) == "N"] <- "pathway.N"
colnames(genome.annotations)[colnames(genome.annotations) == "genome.N"] <- "genome.pathway.N"
colnames(genome.annotations)[colnames(genome.annotations) == "perc.complete"] <- "pathway.perc.complete"

genome.annotations <- merge(x = genome.annotations, y = module.ko.key, by = "ko", all.x = T, all.y = F, allow.cartesian = TRUE)
genome.annotations <- merge(x = genome.annotations, y = genome.modules, by = "module", all.x = T, all.y = F)
colnames(genome.annotations)[colnames(genome.annotations) == "N"] <- "module.N"
colnames(genome.annotations)[colnames(genome.annotations) == "genome.N"] <- "genome.module.N"
colnames(genome.annotations)[colnames(genome.annotations) == "perc.complete"] <- "module.perc.complete"

# this is giving a segfault error on the server for ONE file: ME2012-05-05_3300042888_group4_bin209.kofamscan.tsv.gz
# genome.annotations <- genome.annotations[ ,.(gene, threshold, bit.score, e.value, signif, relaxed.signif.6, 
#                                              ko, ko.description, 
#                                              module, module.description, module.N, genome.module.N, module.perc.complete,
#                                              path, pathway.description, pathway.N, genome.pathway.N, pathway.perc.complete)]
# replacing with base R solves it. left a github issue, I think it's a data.table compilation problem showing up with the updated ubuntu servers (works on my mac, too)
genome.annotations <- genome.annotations[ ,c("gene", "threshold", "bit.score", "e.value", "signif", "relaxed.signif.6", 
                                             "ko", "ko.description", 
                                             "module", "module.description", "module.N", "genome.module.N", "module.perc.complete",
                                             "path", "pathway.description", "pathway.N", "genome.pathway.N", "pathway.perc.complete")]

# ---- one line per gene, all pathways and modules----
# given the variability in pathway granularities, hard to choose a completeness cutoff or "best" pathway criteria
# for now just list all the gene's options

# grep(pattern = ";", x = genome.pathways$pathway.description) # there are no semicolons in the pathway names
# grep(pattern = ";", x = genome.modules$module.description) # or the module names

genome.annotations.combined <- genome.annotations[ , .(path = paste(unique(path), collapse = "; "),
                                                       pathway.description = paste(unique(paste0(pathway.description, " (", genome.pathway.N, "/",pathway.N,"=",round(pathway.perc.complete,0),"%)")), collapse = "; "),
                                                       module = paste(unique(module), collapse = "; "),
                                                       module.description = paste(unique(paste0(module.description, " (", genome.module.N, "/",module.N,"=",round(module.perc.complete,0),"%)")), collapse = "; ")), 
                                                   by = .(gene, threshold, bit.score, e.value, signif, relaxed.signif.6, ko, ko.description)]
genome.annotations.combined[module.description == "NA (NA/NA=NA%)", module.description := paste0("(",ko,": ",ko.description,")")]
genome.annotations.combined[pathway.description == "NA (NA/NA=NA%)", pathway.description := paste0("(",ko,": ",ko.description,")")]

# ---- export data ----

fwrite(x = genome.annotations, file = output.annotations.file.one.ann.per.line, sep = "\t", nThread = num.threads, compress = "gzip")

fwrite(x = genome.annotations.combined, file = output.annotations.file.one.gene.per.line, sep = "\t", nThread = num.threads, compress = "gzip")
