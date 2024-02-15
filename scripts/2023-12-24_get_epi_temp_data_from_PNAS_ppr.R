# RRR 
# data for fig 4b
# just cite the PNAS paper and pull it's data, because it combined many datasets and calculated the epilimnion depth already

library(data.table)
library(lubridate)

temp <- readRDS(file = "pop/data/environmental_data/Robin-Refined/water-temp/10-water_temp_profiles_with_thermocline_stats.rds")
temp <- as.data.table(temp)
colnames(temp)

# require at least 5 epi and 5 hypo measurements to include the profile
index <- grep(pattern = "[[:digit:]]m$", x = colnames(temp), value = F)
index.epi <- index[1:20]
index.hypo <- index[21:52]
epi.mat <- temp[ ,..index.epi]
epi.mat <- as.matrix(epi.mat)
epi.mat <- !is.na(epi.mat)
summary(rowSums(epi.mat, na.rm = T))
hypo.mat <- temp[ ,..index.hypo]
hypo.mat <- as.matrix(hypo.mat)
hypo.mat <- !is.na(hypo.mat)
summary(rowSums(hypo.mat, na.rm = T))
index.keep.epi <- rowSums(epi.mat, na.rm = T) >= 5
index.keep.hypo <- rowSums(hypo.mat, na.rm = T) >= 5
temp <- temp[index.keep.epi == TRUE & index.keep.hypo == TRUE]

# get into long format
temp <- melt(data = temp, id.vars = c("Sample.DateTime","Year","Month","Day","Hour","Minute","Source.Water.Temp"), measure.vars = c(grep(pattern = "[[:digit:]]m$", x = colnames(temp), value = T),"Epi.Mean.Temp.C","Hypo.Mean.Temp.C"), variable.name = "Depth.m", value.name = "Temp.C")
temp[ ,yDay := yday(Sample.DateTime)]

# get info for the plot:
temp <- temp[Depth.m == "Epi.Mean.Temp.C"]

# export simple csv to pair with the plot
temp <- temp[ ,.(Year,yDay,Temp.C)]

fwrite(x = temp, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/epilimnion_mean_temperatures.csv")
