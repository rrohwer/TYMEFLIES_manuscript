# RRR
# format the zoop data for fig 4 plot
# pull from spiny anoxia project because the biomass calculations were already applied
# plus new zoop data isn't on EDI yet anyway

zoops <- readRDS(file = "~/Desktop/spinyAnoxia/data_processed/0q_zoop_abundances.rds") # this includes Jake's biomass calculations
zoops <- as.data.table(zoops)

# remove large predatory zoops that also are unreliably caught in the net
unique(zoops$species_name)
zoops <- zoops[species_name != "Bythotrephes Longimanus" & species_name != "Leptodora Kindti"]
zoops[is.na(Biomass.mg.L), Biomass.mg.L := 0]

tot.zoop <- zoops[ , .(Biomass.mg.L = mean(Biomass.mg.L)), by = .(year4,yday)]

year.means <- tot.zoop[ ,.(Mean.Biomass.mg.L = mean(Biomass.mg.L)), by = .(year4)]

fwrite(x = year.means, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/zoop_mean_annual_biomass.csv")