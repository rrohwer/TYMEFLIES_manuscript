# RRR
# What's the overlap btwn seasonal and long-term trends?

# OK these plots don't really work, I didn't think through it fully
# and the one that does work is kind of stupid, takes up a lot of space to show something simple
# will try sideways bars instead. 


library(data.table)

seas <- fread(file = "data/2023-12-06_abundance_seasonality_analysis/genome_fft_stats_and_bloom_diversity.tsv.gz")
lt <- fread(file = "data/2023-12-04_long-term_change_classifier/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")

colnames(seas)
unique(seas$blooms.are)

seas[ ,`:=`(div.not.seasonal = FALSE)]
seas[nuc.div.is.seasonal == FALSE, `:=`(div.not.seasonal = TRUE)]
seas[ ,`:=`(abund.not.seasonal = FALSE)]
seas[abund.is.seasonal == FALSE, `:=`(abund.not.seasonal = TRUE)]

seas <- seas[ ,.( abund.is.seasonal, nuc.div.is.seasonal, div.not.seasonal, abund.not.seasonal)]

fit <- euler(seas, shape = "circle")
plot(fit, fills = c("#e79493","#a9adff","grey80","grey80"))


# Better to use all the same kind of data in the same plot
# Because saying "we see these lt patterns overlayed over the seasonal ones"
# So best to use the seasonal patterns id'ed from the same SNV data
lt[ ,is.gradual := FALSE]
lt[ Classified.LT.Change == "gradual", is.gradual := TRUE]

lt[ ,is.step := FALSE]
lt[ Classified.LT.Change == "step", is.step := TRUE]

lt[ ,is.dist := FALSE]
lt[ Classified.LT.Change == "disturbance", is.dist := TRUE]

lt[ ,no.LT := FALSE]
lt[ Classified.LT.Change == "none", no.LT := TRUE]

lt[ ,total := TRUE]

lt <- lt[ ,.(Classified.Seasonal, is.gradual, is.step, is.dist, no.LT,total)]

fit <- euler(lt, shape = "circle")
plot(fit, quantities = T)

lt <- lt[ ,.(Classified.Seasonal, is.gradual, is.step, is.dist,total)]

fit <- euler(lt, shape = "circle")
plot(fit, quantities = T, fills = c("grey80","cadetblue2","magenta3","orange2", adjustcolor("white",0)), labels = c(label = c("Seasonal","Gradual Change","Disturbance/Resilience","Step Change", "Total Genomes")))





lt <- fread(file = "data/2023-12-04_long-term_change_classifier/final_files/all_genomes_time_decay_stats-all_SNVs.tsv.gz")

lt[ ,no.LT := FALSE]
lt[ Classified.LT.Change == "none", no.LT := TRUE]
lt[ ,yes.LT := FALSE]
lt[ Classified.LT.Change != "none", no.LT := TRUE]
lt[ ,no.seasonal := FALSE]
lt[ Classified.Seasonal == FALSE, no.seasonal := TRUE]

lt <- lt[ ,.(no.LT, yes.LT, no.seasonal, Classified.Seasonal)]
lt <- as.data.frame(lt)

fit <- euler(lt, shape = "ellipse")
plot(fit)
