# RRR
# get biomass data for phytoplankton
# this time download from EDI
# do have to convert it to biomass though, but check previous script
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.88.31

library(data.table)
library(lubridate)

# ---- download from EDI ----

# Package ID: knb-lter-ntl.88.31 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Phytoplankton - Madison Lakes Area 1995 - current.
# Data set creator:  John Magnuson - University of Wisconsin 
# Data set creator:  Stephen Carpenter - University of Wisconsin 
# Data set creator:  Emily Stanley - University of Wisconsin 
# Contact:    -  NTL LTER  - ntl.infomgr@gmail.com
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/88/31/f2de15b2fff6ae962a04c150c0a1c510" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lakeid",     
                 "year4",     
                 "sampledate",     
                 "sta",     
                 "depth_range",     
                 "division",     
                 "taxa_name",     
                 "gald",     
                 "cells_per_nu",     
                 "nu_per_ml",     
                 "cells_per_ml",     
                 "biovolume_conc",     
                 "biomass_conc",     
                 "relative_total_biovolume",     
                 "genus"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$lakeid)!="factor") dt1$lakeid<- as.factor(dt1$lakeid)
if (class(dt1$year4)=="factor") dt1$year4 <-as.numeric(levels(dt1$year4))[as.integer(dt1$year4) ]               
if (class(dt1$year4)=="character") dt1$year4 <-as.numeric(dt1$year4)                                   
# attempting to convert dt1$sampledate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1sampledate<-as.Date(dt1$sampledate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1sampledate) == length(tmp1sampledate[!is.na(tmp1sampledate)])){dt1$sampledate <- tmp1sampledate } else {print("Date conversion failed for dt1$sampledate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1sampledate) 
if (class(dt1$sta)!="factor") dt1$sta<- as.factor(dt1$sta)
if (class(dt1$depth_range)!="factor") dt1$depth_range<- as.factor(dt1$depth_range)
if (class(dt1$division)!="factor") dt1$division<- as.factor(dt1$division)
if (class(dt1$taxa_name)!="factor") dt1$taxa_name<- as.factor(dt1$taxa_name)
if (class(dt1$gald)=="factor") dt1$gald <-as.numeric(levels(dt1$gald))[as.integer(dt1$gald) ]               
if (class(dt1$gald)=="character") dt1$gald <-as.numeric(dt1$gald)
if (class(dt1$cells_per_nu)=="factor") dt1$cells_per_nu <-as.numeric(levels(dt1$cells_per_nu))[as.integer(dt1$cells_per_nu) ]               
if (class(dt1$cells_per_nu)=="character") dt1$cells_per_nu <-as.numeric(dt1$cells_per_nu)
if (class(dt1$nu_per_ml)=="factor") dt1$nu_per_ml <-as.numeric(levels(dt1$nu_per_ml))[as.integer(dt1$nu_per_ml) ]               
if (class(dt1$nu_per_ml)=="character") dt1$nu_per_ml <-as.numeric(dt1$nu_per_ml)
if (class(dt1$cells_per_ml)=="factor") dt1$cells_per_ml <-as.numeric(levels(dt1$cells_per_ml))[as.integer(dt1$cells_per_ml) ]               
if (class(dt1$cells_per_ml)=="character") dt1$cells_per_ml <-as.numeric(dt1$cells_per_ml)
if (class(dt1$biovolume_conc)=="factor") dt1$biovolume_conc <-as.numeric(levels(dt1$biovolume_conc))[as.integer(dt1$biovolume_conc) ]               
if (class(dt1$biovolume_conc)=="character") dt1$biovolume_conc <-as.numeric(dt1$biovolume_conc)
if (class(dt1$biomass_conc)=="factor") dt1$biomass_conc <-as.numeric(levels(dt1$biomass_conc))[as.integer(dt1$biomass_conc) ]               
if (class(dt1$biomass_conc)=="character") dt1$biomass_conc <-as.numeric(dt1$biomass_conc)
if (class(dt1$relative_total_biovolume)=="factor") dt1$relative_total_biovolume <-as.numeric(levels(dt1$relative_total_biovolume))[as.integer(dt1$relative_total_biovolume) ]               
if (class(dt1$relative_total_biovolume)=="character") dt1$relative_total_biovolume <-as.numeric(dt1$relative_total_biovolume)
if (class(dt1$genus)!="factor") dt1$genus<- as.factor(dt1$genus)


# ---- format ----

phyto <- as.data.table(dt1)

phyto <- phyto[lakeid == "ME", ]
head(phyto)

unique(phyto$sta) # only 1
unique(phyto$depth_range) # 0-8 and 0-2 m pulls
colnames(phyto)

phyto <- phyto[ ,-c(1,4)]
phyto$yday <- yday(phyto$sampledate)

# what are these columns??
# nu = natural unit, (i.e., colonies, filaments, or single cells)
# note: Multiple entries for the same species on the same date may be due to different variants or vegetative states - (e.g., colonial or attached vs. free cell.)
# note: one million cubic Micrometers of biovolume Per mL of water are equal to a biovolume concentration of one mm^3/mL == 1 mg/L (if cell density ~ water)
# A tube sampler is used for the 0-8 m Lake Mendota samples
# gald = avg greatest axial linear dimension, um
# cells_per_nu =  per natural unit 
# biovolume concentration = micrometerCubedPerMilliliter
# biomass conc = milligramsPerLiter
# rel total biovolume = perc

colnames(phyto)[colnames(phyto) == "depth_range"] <- "Sample.Depth.m"
colnames(phyto)[colnames(phyto) == "gald"] <- "Size.length.um"
colnames(phyto)[colnames(phyto) == "cells_per_nu"] <- "Size.cells.unit"
colnames(phyto)[colnames(phyto) == "nu_per_ml"] <- "Count.units.mL"
colnames(phyto)[colnames(phyto) == "cells_per_ml"] <- "Count.cells.mL"
colnames(phyto)[colnames(phyto) == "biovolume_conc"] <- "Biovolume.um3.mL"
colnames(phyto)[colnames(phyto) == "biomass_conc"] <- "Biomass.mg.L"

colnames(phyto)

# remove trailing white spaces from names:
grep(pattern = "Microcystis", x = phyto$taxa_name, ignore.case = T, value = T) # e.g.

phyto$taxa_name <- sub(pattern = " $", "", phyto$taxa_name)
phyto$genus <- sub(pattern = " $", "", phyto$genus)
phyto$division <- sub(pattern = " $", "", phyto$division)

phyto$taxa_name <- sub(pattern = " $", "", phyto$taxa_name)
phyto$genus <- sub(pattern = " $", "", phyto$genus)
phyto$division <- sub(pattern = " $", "", phyto$division)

phyto$taxa_name <- sub(pattern = " $", "", phyto$taxa_name)
phyto$genus <- sub(pattern = " $", "", phyto$genus)
phyto$division <- sub(pattern = " $", "", phyto$division)

# ---- get annual averages ----

# remove outlandish values
phyto <- phyto[Biomass.mg.L < 80]

# first get sum on each day
daily <- phyto[ , .(Biomass.mg.L = sum(Biomass.mg.L)), by = .(year4, yday)] 

# check it
plot(x = daily$yday, y = daily$Biomass.mg.L, type = "n")
for (yr in unique(daily$year4)){
  points(x = daily[year4 == yr, yday], y = daily[year4 == yr, Biomass.mg.L], pch = 16, col = adjustcolor("grey40",.5))
  lines(x = daily[year4 == yr, yday], y = daily[year4 == yr, Biomass.mg.L], col = adjustcolor("grey40",.5))
}
points(x = daily[year4 == 2012, yday], y = daily[year4 == 2012, Biomass.mg.L], pch = 16, col = adjustcolor("red2",.5))
lines(x = daily[year4 == 2012, yday], y = daily[year4 == 2012, Biomass.mg.L], col = adjustcolor("red2",.5), lwd = 3)

# now get annual means
yearly <- daily[ , .(Mean.Biomass.mg.L = mean(Biomass.mg.L, na.rm = T)), by = .(year4)]

# check it
boxplot(yearly$Mean.Biomass.mg.L)
points(x = 1, y = yearly[year4 == 2012, Mean.Biomass.mg.L], col = "red2", pch = 16)

# ---- export data ----

fwrite(x = yearly, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/phyto_mean_biomass.csv")

