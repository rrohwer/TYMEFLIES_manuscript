# RRR
# Get which pathways are in the genome.
# Decisions:
# 1. What is a significant KO annotation?
#       If use the KEGG star (default kofamscan significance) there are basically no complete pathways in the acI
#       Katy and Pedro use e-value as a secondary cutoff (1e-5)
#       Anvi'o relaxes the threshold for bit scores to include. https://anvio.org/help/main/programs/anvi-run-kegg-kofams/
#       Looking at the defaults, the threshold is a manual decision for each hmmer family about how good it should match (how good a bit score)
#       And the default e-value is high, just < .01
#       So defaults are relying more on threshold/bit score than e-value to say it's significant. 
#       So I will do that too. keep e-value cutoff at < .01 and give second-tier classifications by relaxing threshold.
# 2. What is a complete pathway?
#       Partly depends on how many enzymes are in the pathway!
# 3. Which pathways reference to use? 
#       For now, just using the default kegg-defined pathways. 
#       In future, could sub in other people's specific pathways instead, like the keggdecoder ones.

# ---- set-up ----

library(data.table)
library(ggplot2)
library(lubridate)

pathways.key <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/ko_pathway_links.txt", header = F)

module.ko.key <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/ko_module_links.txt", header = F)

module.path.key <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/module_pathway_links.txt", header = F)

pathways.names <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/pathway_list.txt", header = F, sep = "\t")

ko.names <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/ko_list.txt", header = F, sep = "\t")

module.names <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/module_list.txt", header = F, sep = "\t")

genome.annotations <- fread(file = "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.kofamscan.tsv")

plot.folder <- "figures/2023-09-25_exploring_KOs_and_pathways/"

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

# ---- which KO's are in the genome ----

genome.annotations$default.signif <- FALSE
genome.annotations[signif == "*", default.signif := TRUE]

genome.annotations$relaxed.signif.9 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .9), relaxed.signif.9 := TRUE]

genome.annotations$relaxed.signif.8 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .8), relaxed.signif.8 := TRUE]

genome.annotations$relaxed.signif.7 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .7), relaxed.signif.7 := TRUE]

genome.annotations$relaxed.signif.6 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .6), relaxed.signif.6 := TRUE]

genome.annotations$relaxed.signif.5 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .5), relaxed.signif.5 := TRUE]

genome.annotations$relaxed.signif.4 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .4), relaxed.signif.4 := TRUE]

genome.annotations$relaxed.signif.3 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .3), relaxed.signif.3 := TRUE]

genome.annotations$relaxed.signif.2 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .2), relaxed.signif.2 := TRUE]

genome.annotations$relaxed.signif.1 <- FALSE
genome.annotations[e.value < 0.01 & bit.score > (threshold * .1), relaxed.signif.1 := TRUE]

genome.annotations[ ,signif := NULL]

# Remove the massive pathways ----

selected.pathways <- pathways.key[ , .N, by = .(path)]
selected.pathways <- merge(x = pathways.names, y = selected.pathways, by = "path")
selected.pathways <- selected.pathways[order(N, decreasing = T)]
# View(selected.pathways)
selected.pathways <- selected.pathways[-(1:5)]

pathways.key <- merge(x = selected.pathways, y = pathways.key, by = "path", all = F)

# no massive modules
selected.modules <- module.ko.key[ , .N, by = .(module)]
selected.modules[order(N, decreasing = T)]

module.ko.key <- merge(x = module.ko.key, y = module.names, by = "module")

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

genome.pathways <- pathways.key[ , .N, by = .(path, pathway.description)]
genome.pathways$genome.N.default <- 0
genome.pathways$genome.N.9 <- 0
genome.pathways$genome.N.8 <- 0
genome.pathways$genome.N.7 <- 0
genome.pathways$genome.N.6 <- 0
genome.pathways$genome.N.5 <- 0
genome.pathways$genome.N.4 <- 0
genome.pathways$genome.N.3 <- 0
genome.pathways$genome.N.2 <- 0
genome.pathways$genome.N.1 <- 0

for (p in genome.pathways$path){
  path.kos <- pathways.list[[p]]
  
  x <- intersect(x = path.kos, genome.annotations[default.signif == TRUE , ko])
  genome.pathways[path == p, genome.N.default := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.9 == TRUE , ko])
  genome.pathways[path == p, genome.N.9 := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.8 == TRUE , ko])
  genome.pathways[path == p, genome.N.8 := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.7 == TRUE , ko])
  genome.pathways[path == p, genome.N.7 := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.6 == TRUE , ko])
  genome.pathways[path == p, genome.N.6 := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.5 == TRUE , ko])
  genome.pathways[path == p, genome.N.5 := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.4 == TRUE , ko])
  genome.pathways[path == p, genome.N.4 := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.3 == TRUE , ko])
  genome.pathways[path == p, genome.N.3 := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.2 == TRUE , ko])
  genome.pathways[path == p, genome.N.2 := length(x)]
  
  x <- intersect(x = path.kos, genome.annotations[relaxed.signif.1 == TRUE , ko])
  genome.pathways[path == p, genome.N.1 := length(x)]
}

genome.pathways[ ,perc.complete.default := genome.N.default / N * 100]
genome.pathways[ ,perc.complete.9 := genome.N.9 / N * 100]
genome.pathways[ ,perc.complete.8 := genome.N.8 / N * 100]
genome.pathways[ ,perc.complete.7 := genome.N.7 / N * 100]
genome.pathways[ ,perc.complete.6 := genome.N.6 / N * 100]
genome.pathways[ ,perc.complete.5 := genome.N.5 / N * 100]
genome.pathways[ ,perc.complete.4 := genome.N.4 / N * 100]
genome.pathways[ ,perc.complete.3 := genome.N.3 / N * 100]
genome.pathways[ ,perc.complete.2 := genome.N.2 / N * 100]
genome.pathways[ ,perc.complete.1 := genome.N.1 / N * 100]

summary(genome.pathways)

# ---- which modules are in the genome? ----
# for each pathway, list its KOs 

modules.list <- list(NULL)
for (m in unique(module.ko.key$module)){
  k <- module.ko.key[module == m]
  k <- list(k$ko)
  names(k) <- m
  modules.list <- c(modules.list, k)
}

# for each module, how many of its KOs are present?

genome.modules <- module.ko.key[ , .N, by = .(module, module.description)]
genome.modules$genome.N.default <- 0
genome.modules$genome.N.9 <- 0
genome.modules$genome.N.8 <- 0
genome.modules$genome.N.7 <- 0
genome.modules$genome.N.6 <- 0
genome.modules$genome.N.5 <- 0
genome.modules$genome.N.4 <- 0
genome.modules$genome.N.3 <- 0
genome.modules$genome.N.2 <- 0
genome.modules$genome.N.1 <- 0

for (m in genome.modules$module){
  module.kos <- modules.list[[m]]
  
  x <- intersect(x = module.kos, genome.annotations[default.signif == TRUE , ko])
  genome.modules[module == m, genome.N.default := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.9 == TRUE , ko])
  genome.modules[module == m, genome.N.9 := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.8 == TRUE , ko])
  genome.modules[module == m, genome.N.8 := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.7 == TRUE , ko])
  genome.modules[module == m, genome.N.7 := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.6 == TRUE , ko])
  genome.modules[module == m, genome.N.6 := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.5 == TRUE , ko])
  genome.modules[module == m, genome.N.5 := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.4 == TRUE , ko])
  genome.modules[module == m, genome.N.4 := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.3 == TRUE , ko])
  genome.modules[module == m, genome.N.3 := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.2 == TRUE , ko])
  genome.modules[module == m, genome.N.2 := length(x)]
  
  x <- intersect(x = module.kos, genome.annotations[relaxed.signif.1 == TRUE , ko])
  genome.modules[module == m, genome.N.1 := length(x)]
}

genome.modules[ ,perc.complete.default := genome.N.default / N * 100]
genome.modules[ ,perc.complete.9 := genome.N.9 / N * 100]
genome.modules[ ,perc.complete.8 := genome.N.8 / N * 100]
genome.modules[ ,perc.complete.7 := genome.N.7 / N * 100]
genome.modules[ ,perc.complete.6 := genome.N.6 / N * 100]
genome.modules[ ,perc.complete.5 := genome.N.5 / N * 100]
genome.modules[ ,perc.complete.4 := genome.N.4 / N * 100]
genome.modules[ ,perc.complete.3 := genome.N.3 / N * 100]
genome.modules[ ,perc.complete.2 := genome.N.2 / N * 100]
genome.modules[ ,perc.complete.1 := genome.N.1 / N * 100]

summary(genome.modules)

# How does pathway completeness increase with relaxed thresholds? ----

genome.pathways.long <- melt(data = genome.pathways, id.vars = c("path", "pathway.description", "N"), 
                             measure.vars = c("perc.complete.default","perc.complete.9","perc.complete.8","perc.complete.7","perc.complete.6","perc.complete.5","perc.complete.4","perc.complete.3","perc.complete.2","perc.complete.1"), 
                             variable.name = "Threshold.perc", value.name = "perc.completeness")
genome.pathways.long[ ,Threshold.perc := sub(pattern = "default",replacement = "10",x = Threshold.perc)]
genome.pathways.long[ ,Threshold.perc := sub(pattern = "perc\\.complete\\.",replacement = "",x = Threshold.perc)]
genome.pathways.long[ ,Threshold.perc := as.numeric(Threshold.perc) * 10]

# How does module completeness increase with relaxed thresholds? ----

genome.modules.long <- melt(data = genome.modules, id.vars = c("module", "module.description", "N"), 
                             measure.vars = c("perc.complete.default","perc.complete.9","perc.complete.8","perc.complete.7","perc.complete.6","perc.complete.5","perc.complete.4","perc.complete.3","perc.complete.2","perc.complete.1"), 
                             variable.name = "Threshold.perc", value.name = "perc.completeness")
genome.modules.long[ ,Threshold.perc := sub(pattern = "default",replacement = "10",x = Threshold.perc)]
genome.modules.long[ ,Threshold.perc := sub(pattern = "perc\\.complete\\.",replacement = "",x = Threshold.perc)]
genome.modules.long[ ,Threshold.perc := as.numeric(Threshold.perc) * 10]

# really should subset to how much the existing modules are completed at relaxed thresholds...

genome.modules.subset <- genome.modules[genome.N.default > 0]

genome.modules.long.subset <- melt(data = genome.modules.subset, id.vars = c("module", "module.description", "N"), 
                            measure.vars = c("perc.complete.default","perc.complete.9","perc.complete.8","perc.complete.7","perc.complete.6","perc.complete.5","perc.complete.4","perc.complete.3","perc.complete.2","perc.complete.1"), 
                            variable.name = "Threshold.perc", value.name = "perc.completeness")
genome.modules.long.subset[ ,Threshold.perc := sub(pattern = "default",replacement = "10",x = Threshold.perc)]
genome.modules.long.subset[ ,Threshold.perc := sub(pattern = "perc\\.complete\\.",replacement = "",x = Threshold.perc)]
genome.modules.long.subset[ ,Threshold.perc := as.numeric(Threshold.perc) * 10]


# pathway plots ----

pdf(file = file.path(plot.folder,"KO_annotation_threshold_choice.pdf"), width = 9, height = 7)
ggplot(data = genome.pathways.long, aes(y = perc.completeness, x = Threshold.perc))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_boxplot(aes(group = Threshold.perc), outlier.shape = NA)+
  geom_jitter(aes(color = N))+
  scale_color_gradient(low = adjustcolor("grey",.2), high = "magenta", name = "Number\nof KOs in\nPathway")+
  labs(title = "How does pathway completeness increase with relaxed annotation threshold?", 
       subtitle = "Each dot is a KEGG pathway.")+
  ylab(label = "Pathway Completeness (% total KOs in pathway)")+
  scale_x_continuous(name = "Percent of Default Threshold Applied (%)", breaks = seq.int(from = 100, to = 10, by = -10))
dev.off()

pdf(file = file.path(plot.folder,"KO_annotation_threshold_choice-subset_smaller_pathways.pdf"), width = 9, height = 7)
ggplot(data = genome.pathways.long[N < 30], aes(y = perc.completeness, x = Threshold.perc))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_boxplot(aes(group = Threshold.perc), outlier.shape = NA)+
  geom_jitter(aes(color = N))+
  scale_color_gradient(low = adjustcolor("grey",.2), high = "magenta", name = "Number\nof KOs in\nPathway")+
  labs(title = "How does pathway completeness increase with relaxed annotation threshold?", 
       subtitle = "Each dot is a KEGG pathway.")+
  ylab(label = "Pathway Completeness (% total KOs in pathway)")+
  scale_x_continuous(name = "Percent of Default Threshold Applied (%)", breaks = seq.int(from = 100, to = 10, by = -10))
dev.off()

pdf(file = file.path(plot.folder,"KO_annotation_threshold_choice-no_zeros.pdf"), width = 9, height = 7)
ggplot(data = genome.pathways.long[perc.completeness > 0], aes(y = perc.completeness, x = Threshold.perc))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_boxplot(aes(group = Threshold.perc), outlier.shape = NA)+
  geom_jitter(aes(color = N))+
  scale_color_gradient(low = adjustcolor("grey",.2), high = "magenta", name = "Number\nof KOs in\nPathway")+
  labs(title = "How does pathway completeness increase with relaxed annotation threshold?", 
       subtitle = "Each dot is a KEGG pathway.")+
  ylab(label = "Pathway Completeness (% total KOs in pathway)")+
  scale_x_continuous(name = "Percent of Default Threshold Applied (%)", breaks = seq.int(from = 100, to = 10, by = -10))
dev.off()

# honestly I'd probably choose .6 as the cutoff, not .5

pdf(file = file.path(plot.folder,"Pathway_completeness_in_acI-B.pdf"), width = 8, height = 20)
ggplot(data = genome.pathways.long[Threshold.perc == 60 & perc.completeness > 5], aes(x = perc.completeness, y = pathway.description))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank())+
  geom_col(aes(fill = perc.completeness), col = "thistle4")+
  scale_fill_gradient(low = adjustcolor("thistle3",.2), high = "thistle3", guide = F)+
  scale_x_continuous(expand = c(0,0))+
  labs(title = "The 136 pathways > 5% complete in acI-B",
       subtitle = "not obvious how to choose a completeness cutoff...")+
  xlab(label = "Pathway Completeness (%)")
dev.off()

# These pathways are not looking very complete... but that's because they're pretty broad
# right, these are summary categories, not distinct pathways here

# module plots ----

pdf(file = file.path(plot.folder,"KO_annotation_threshold_choice-modules.pdf"), width = 9, height = 7)
ggplot(data = genome.modules.long, aes(y = perc.completeness, x = Threshold.perc))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_boxplot(aes(group = Threshold.perc), outlier.shape = NA)+
  geom_jitter(aes(color = N))+
  scale_color_gradient(low = adjustcolor("grey",.2), high = "magenta", name = "Number\nof KOs in\nModule")+
  labs(title = "How does pathway completeness increase with relaxed annotation threshold?", 
       subtitle = "Each dot is a KEGG module")+
  ylab(label = "Module Completeness (% total KOs in module)")+
  scale_x_continuous(name = "Percent of Default Threshold Applied (%)", breaks = seq.int(from = 100, to = 10, by = -10))
dev.off()

pdf(file = file.path(plot.folder,"KO_annotation_threshold_choice-no_zeros-modules.pdf"), width = 9, height = 7)
ggplot(data = genome.modules.long[perc.completeness > 0], aes(y = perc.completeness, x = Threshold.perc))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_boxplot(aes(group = Threshold.perc), outlier.shape = NA)+
  geom_jitter(aes(color = N))+
  scale_color_gradient(low = adjustcolor("grey",.2), high = "magenta", name = "Number\nof KOs in\nModule")+
  labs(title = "How does module completeness increase with relaxed annotation threshold?", 
       subtitle = "Each dot is a KEGG module")+
  ylab(label = "Module Completeness (% total KOs in module)")+
  scale_x_continuous(name = "Percent of Default Threshold Applied (%)", breaks = seq.int(from = 100, to = 10, by = -10))
dev.off()

# look at how existing modules grow, not completeness of all modules

pdf(file = file.path(plot.folder,"KO_annotation_threshold_choice-modules-that-exist-at-default-cutoff-only.pdf"), width = 9, height = 7)
ggplot(data = genome.modules.long.subset, aes(y = perc.completeness, x = Threshold.perc))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_boxplot(aes(group = Threshold.perc), outlier.shape = NA)+
  geom_jitter(aes(color = N))+
  scale_color_gradient(low = adjustcolor("grey",.2), high = "magenta", name = "Number\nof KOs in\nModule")+
  labs(title = "How does pathway completeness increase with relaxed annotation threshold?", 
       subtitle = "Each dot is a KEGG module")+
  ylab(label = "Module Completeness (% total KOs in module)")+
  scale_x_continuous(name = "Percent of Default Threshold Applied (%)", breaks = seq.int(from = 100, to = 10, by = -10))
dev.off()

# cutoff choice still not obvious to me. as reminder, anvi'o uses .5


pdf(file = file.path(plot.folder,"Module_completeness_in_acI-B.pdf"), width = 8, height = 20)
ggplot(data = genome.modules.long[Threshold.perc == 60 & perc.completeness > 5], aes(x = perc.completeness, y = module.description))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank())+
  geom_col(aes(fill = perc.completeness), col = "thistle4")+
  scale_fill_gradient(low = adjustcolor("thistle3",.2), high = "thistle3", guide = F)+
  scale_x_continuous(expand = c(0,0))+
  labs(title = "The 201 modules > 5% complete in acI-B",
       subtitle = "not obvious how to choose a completeness cutoff...")+
  xlab(label = "Module Completeness (%)")
dev.off()









# ---- which pathways do the genome's KOs belong to ----

# only include genes with known annotation and known pathways

genome.annotations <- merge(x = genome.annotations, y = pathways.key, by = "ko", all.x = F, all.y = F)

# only include pathways that are at least 5% complete

genome.annotations <- merge(x = genome.annotations, y = genome.pathways[perc.complete.6 > 5], by = c("path","pathway.description","N"), all.x = F, all.y = F)


# ---- pull in genes under selection data ----



genes <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz")
genes <- genes[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman)]

genes <- merge(x = genome.annotations, y = genes, by = "gene", all.x = F, all.y = T, allow.cartesian = T)

colnames(genes)[colnames(genes) == "N"] <- "N.pathway"
genes <- genes[ , .N, by = .(date, pathway.description, N.pathway)]
genes <- genes[order(N.pathway, decreasing = T)]

# genes <- dcast(data = genes, formula = pathway.description ~ date, value.var = "N")
# genes <- as.matrix(genes, rownames = "pathway.description")

# my.colors <- colorRampPalette(colors = c("white","red4"))

genes[ ,date := parse_date_time(date, orders = "ymd")]

pdf(file = file.path(plot.folder,"genes_under_selection_example-acI-B.pdf"), width = 24, height = 24)
ggplot(data = genes, aes(x = date, y = pathway.description, color = N))+
  geom_point(aes(size = N), alpha = .5)+
  theme_bw()+
  scale_color_gradient(low = "grey90", high = "red4")
dev.off()


