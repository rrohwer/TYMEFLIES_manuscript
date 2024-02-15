# RRR
# What's the overlap btwn seasonal and long-term trends?

library(data.table)

lt <- fread(file = "data/2023-12-04_long-term_change_classifier/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")

lt[Classified.Seasonal == TRUE, temp.col := "seasonal"]
lt[Classified.Seasonal == FALSE, temp.col := "not.seasonal"]
lt <- dcast(data = lt, formula = Classified.LT.Change ~ temp.col, fun.aggregate = length)

lt <- as.matrix(lt, rownames = T)
lt <- lt[c(1,4,2),c(2,1)]
lt <- t(lt)

# # check
# barplot(lt, horiz = T, legend.text = T)

saveRDS(object = lt, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_3_data/invisible_present_barplot_matrix.rds")
