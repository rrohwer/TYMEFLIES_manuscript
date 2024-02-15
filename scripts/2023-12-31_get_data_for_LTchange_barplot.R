# RRR
# What's the spread of LT change over different taxa?
# plot by taxon mean abundances

library(data.table)

lt <- fread(file = "data/2023-12-04_long-term_change_classifier/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")
tax <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_genome.tsv")

# ---- format ----

tax <- tax[ ,.(genome, domain, phylum, class, order, family, genus, species, mean.abund)]

lt <- lt[Classified.LT.Change != "none" ,.(genome,Classified.LT.Change)]

lt <- merge(x = lt, y = tax, by = "genome", all.x = T, all.y = F)

lt[ ,Classified.LT.Change := factor(Classified.LT.Change, levels = c("disturbance","step","gradual"), ordered = T)]
lt <- lt[order(phylum,Classified.LT.Change), .(genome,phylum,Classified.LT.Change, mean.abund)]
lt <- dcast(data = lt, formula = genome + phylum + Classified.LT.Change ~ phylum, value.var = "mean.abund")
lt <- lt[order(phylum,Classified.LT.Change, Actinobacteriota,Bacteroidota,Chloroflexota,Cyanobacteria,Planctomycetota,Proteobacteria,Verrucomicrobiota)]

by.genome <- as.matrix(lt[ ,-c(2,3)], rownames = T)
by.genome[is.na(by.genome)] <- 0

by.change <- lt[ , .(Actinobacteriota = sum(Actinobacteriota, na.rm = T),
                     Bacteroidota = sum(Bacteroidota, na.rm = T),
                     Chloroflexota = sum(Chloroflexota, na.rm = T),
                     Cyanobacteria = sum(Cyanobacteria, na.rm = T),
                     Planctomycetota = sum(Planctomycetota, na.rm = T),
                     Proteobacteria = sum(Proteobacteria, na.rm = T),
                     Verrucomicrobiota = sum(Verrucomicrobiota, na.rm = T)), by = .(Classified.LT.Change)]
by.change <- as.matrix(by.change, rownames = T)

colindex <- colSums(by.change)
colindex <- order(colindex, decreasing = T)

by.genome <- by.genome[ ,colindex]
by.change <- by.change[ ,colindex]


# # check
# barplot(by.change, col = c("cadetblue2","magenta3","orange2"), beside = F, legend.text = T, las = 2, border = NA)
# barplot(by.genome, col = adjustcolor("white",0), add = T, beside = F, legend.text = F, las = 2)

# ---- save the two tables ----

saveRDS(object = by.genome, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/LT_change_barplot-by_genome_for_outlines.rds")
saveRDS(object = by.change, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/LT_change_barplot-by_change_for_colors.rds")




