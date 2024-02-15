# RRR
# download LTER nutrient data from EDI
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.1.60
# save the TP, SRP, and DOC for my plots 4 f,g,i

library(data.table)
library(lubridate)

# ---- download data from EDI ----

# Package ID: knb-lter-ntl.1.60 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Chemical Limnology of Primary Study Lakes: Nutrients, pH and Carbon 1981 - current.
# Data set creator:  John Magnuson - University of Wisconsin-Madison 
# Data set creator:  Stephen Carpenter - University of Wisconsin-Madison 
# Data set creator:  Emily Stanley - University of Wisconsin-Madison 
# Metadata Provider:  NTL Information Manager - University of Wisconsin-Madison 
# Contact:    -  NTL LTER  - ntl.infomgr@gmail.com
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/1/60/0ff1fd13116d6097376e3745194cdc5f" 
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
                 "daynum",     
                 "sampledate",     
                 "depth",     
                 "rep",     
                 "sta",     
                 "event",     
                 "ph",     
                 "phair",     
                 "alk",     
                 "dic",     
                 "tic",     
                 "doc",     
                 "toc",     
                 "no3no2",     
                 "no2",     
                 "nh4",     
                 "totnf",     
                 "totnuf",     
                 "totpf",     
                 "totpuf",     
                 "drsif",     
                 "brsif",     
                 "brsiuf",     
                 "tpm",     
                 "totnuf_sloh",     
                 "no3no2_sloh",     
                 "nh4_sloh",     
                 "kjdl_n_sloh",     
                 "totpuf_sloh",     
                 "drp_sloh",     
                 "drsif_sloh",     
                 "flagdepth",     
                 "flagph",     
                 "flagphair",     
                 "flagalk",     
                 "flagdic",     
                 "flagtic",     
                 "flagdoc",     
                 "flagtoc",     
                 "flagno3no2",     
                 "flagno2",     
                 "flagnh4",     
                 "flagtotnf",     
                 "flagtotnuf",     
                 "flagtotpf",     
                 "flagtotpuf",     
                 "flagdrsif",     
                 "flagbrsif",     
                 "flagbrsiuf",     
                 "flagtpm",     
                 "flagtotnuf_sloh",     
                 "flagno3no2_sloh",     
                 "flagnh4_sloh",     
                 "flagkjdl_n_sloh",     
                 "flagtotpuf_sloh",     
                 "flagdrp_sloh",     
                 "flagdrsif_sloh"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$lakeid)!="factor") dt1$lakeid<- as.factor(dt1$lakeid)
if (class(dt1$year4)=="factor") dt1$year4 <-as.numeric(levels(dt1$year4))[as.integer(dt1$year4) ]               
if (class(dt1$year4)=="character") dt1$year4 <-as.numeric(dt1$year4)
if (class(dt1$daynum)=="factor") dt1$daynum <-as.numeric(levels(dt1$daynum))[as.integer(dt1$daynum) ]               
if (class(dt1$daynum)=="character") dt1$daynum <-as.numeric(dt1$daynum)                                   
# attempting to convert dt1$sampledate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1sampledate<-as.Date(dt1$sampledate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1sampledate) == length(tmp1sampledate[!is.na(tmp1sampledate)])){dt1$sampledate <- tmp1sampledate } else {print("Date conversion failed for dt1$sampledate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1sampledate) 
if (class(dt1$depth)=="factor") dt1$depth <-as.numeric(levels(dt1$depth))[as.integer(dt1$depth) ]               
if (class(dt1$depth)=="character") dt1$depth <-as.numeric(dt1$depth)
if (class(dt1$rep)!="factor") dt1$rep<- as.factor(dt1$rep)
if (class(dt1$sta)!="factor") dt1$sta<- as.factor(dt1$sta)
if (class(dt1$event)!="factor") dt1$event<- as.factor(dt1$event)
if (class(dt1$ph)=="factor") dt1$ph <-as.numeric(levels(dt1$ph))[as.integer(dt1$ph) ]               
if (class(dt1$ph)=="character") dt1$ph <-as.numeric(dt1$ph)
if (class(dt1$phair)=="factor") dt1$phair <-as.numeric(levels(dt1$phair))[as.integer(dt1$phair) ]               
if (class(dt1$phair)=="character") dt1$phair <-as.numeric(dt1$phair)
if (class(dt1$alk)=="factor") dt1$alk <-as.numeric(levels(dt1$alk))[as.integer(dt1$alk) ]               
if (class(dt1$alk)=="character") dt1$alk <-as.numeric(dt1$alk)
if (class(dt1$dic)=="factor") dt1$dic <-as.numeric(levels(dt1$dic))[as.integer(dt1$dic) ]               
if (class(dt1$dic)=="character") dt1$dic <-as.numeric(dt1$dic)
if (class(dt1$tic)=="factor") dt1$tic <-as.numeric(levels(dt1$tic))[as.integer(dt1$tic) ]               
if (class(dt1$tic)=="character") dt1$tic <-as.numeric(dt1$tic)
if (class(dt1$doc)=="factor") dt1$doc <-as.numeric(levels(dt1$doc))[as.integer(dt1$doc) ]               
if (class(dt1$doc)=="character") dt1$doc <-as.numeric(dt1$doc)
if (class(dt1$toc)=="factor") dt1$toc <-as.numeric(levels(dt1$toc))[as.integer(dt1$toc) ]               
if (class(dt1$toc)=="character") dt1$toc <-as.numeric(dt1$toc)
if (class(dt1$no3no2)=="factor") dt1$no3no2 <-as.numeric(levels(dt1$no3no2))[as.integer(dt1$no3no2) ]               
if (class(dt1$no3no2)=="character") dt1$no3no2 <-as.numeric(dt1$no3no2)
if (class(dt1$no2)=="factor") dt1$no2 <-as.numeric(levels(dt1$no2))[as.integer(dt1$no2) ]               
if (class(dt1$no2)=="character") dt1$no2 <-as.numeric(dt1$no2)
if (class(dt1$nh4)=="factor") dt1$nh4 <-as.numeric(levels(dt1$nh4))[as.integer(dt1$nh4) ]               
if (class(dt1$nh4)=="character") dt1$nh4 <-as.numeric(dt1$nh4)
if (class(dt1$totnf)=="factor") dt1$totnf <-as.numeric(levels(dt1$totnf))[as.integer(dt1$totnf) ]               
if (class(dt1$totnf)=="character") dt1$totnf <-as.numeric(dt1$totnf)
if (class(dt1$totnuf)=="factor") dt1$totnuf <-as.numeric(levels(dt1$totnuf))[as.integer(dt1$totnuf) ]               
if (class(dt1$totnuf)=="character") dt1$totnuf <-as.numeric(dt1$totnuf)
if (class(dt1$totpf)=="factor") dt1$totpf <-as.numeric(levels(dt1$totpf))[as.integer(dt1$totpf) ]               
if (class(dt1$totpf)=="character") dt1$totpf <-as.numeric(dt1$totpf)
if (class(dt1$totpuf)=="factor") dt1$totpuf <-as.numeric(levels(dt1$totpuf))[as.integer(dt1$totpuf) ]               
if (class(dt1$totpuf)=="character") dt1$totpuf <-as.numeric(dt1$totpuf)
if (class(dt1$drsif)=="factor") dt1$drsif <-as.numeric(levels(dt1$drsif))[as.integer(dt1$drsif) ]               
if (class(dt1$drsif)=="character") dt1$drsif <-as.numeric(dt1$drsif)
if (class(dt1$brsif)=="factor") dt1$brsif <-as.numeric(levels(dt1$brsif))[as.integer(dt1$brsif) ]               
if (class(dt1$brsif)=="character") dt1$brsif <-as.numeric(dt1$brsif)
if (class(dt1$brsiuf)=="factor") dt1$brsiuf <-as.numeric(levels(dt1$brsiuf))[as.integer(dt1$brsiuf) ]               
if (class(dt1$brsiuf)=="character") dt1$brsiuf <-as.numeric(dt1$brsiuf)
if (class(dt1$tpm)=="factor") dt1$tpm <-as.numeric(levels(dt1$tpm))[as.integer(dt1$tpm) ]               
if (class(dt1$tpm)=="character") dt1$tpm <-as.numeric(dt1$tpm)
if (class(dt1$totnuf_sloh)=="factor") dt1$totnuf_sloh <-as.numeric(levels(dt1$totnuf_sloh))[as.integer(dt1$totnuf_sloh) ]               
if (class(dt1$totnuf_sloh)=="character") dt1$totnuf_sloh <-as.numeric(dt1$totnuf_sloh)
if (class(dt1$no3no2_sloh)=="factor") dt1$no3no2_sloh <-as.numeric(levels(dt1$no3no2_sloh))[as.integer(dt1$no3no2_sloh) ]               
if (class(dt1$no3no2_sloh)=="character") dt1$no3no2_sloh <-as.numeric(dt1$no3no2_sloh)
if (class(dt1$nh4_sloh)=="factor") dt1$nh4_sloh <-as.numeric(levels(dt1$nh4_sloh))[as.integer(dt1$nh4_sloh) ]               
if (class(dt1$nh4_sloh)=="character") dt1$nh4_sloh <-as.numeric(dt1$nh4_sloh)
if (class(dt1$kjdl_n_sloh)=="factor") dt1$kjdl_n_sloh <-as.numeric(levels(dt1$kjdl_n_sloh))[as.integer(dt1$kjdl_n_sloh) ]               
if (class(dt1$kjdl_n_sloh)=="character") dt1$kjdl_n_sloh <-as.numeric(dt1$kjdl_n_sloh)
if (class(dt1$totpuf_sloh)=="factor") dt1$totpuf_sloh <-as.numeric(levels(dt1$totpuf_sloh))[as.integer(dt1$totpuf_sloh) ]               
if (class(dt1$totpuf_sloh)=="character") dt1$totpuf_sloh <-as.numeric(dt1$totpuf_sloh)
if (class(dt1$drp_sloh)=="factor") dt1$drp_sloh <-as.numeric(levels(dt1$drp_sloh))[as.integer(dt1$drp_sloh) ]               
if (class(dt1$drp_sloh)=="character") dt1$drp_sloh <-as.numeric(dt1$drp_sloh)
if (class(dt1$drsif_sloh)=="factor") dt1$drsif_sloh <-as.numeric(levels(dt1$drsif_sloh))[as.integer(dt1$drsif_sloh) ]               
if (class(dt1$drsif_sloh)=="character") dt1$drsif_sloh <-as.numeric(dt1$drsif_sloh)
if (class(dt1$flagdepth)!="factor") dt1$flagdepth<- as.factor(dt1$flagdepth)
if (class(dt1$flagph)!="factor") dt1$flagph<- as.factor(dt1$flagph)
if (class(dt1$flagphair)!="factor") dt1$flagphair<- as.factor(dt1$flagphair)
if (class(dt1$flagalk)!="factor") dt1$flagalk<- as.factor(dt1$flagalk)
if (class(dt1$flagdic)!="factor") dt1$flagdic<- as.factor(dt1$flagdic)
if (class(dt1$flagtic)!="factor") dt1$flagtic<- as.factor(dt1$flagtic)
if (class(dt1$flagdoc)!="factor") dt1$flagdoc<- as.factor(dt1$flagdoc)
if (class(dt1$flagtoc)!="factor") dt1$flagtoc<- as.factor(dt1$flagtoc)
if (class(dt1$flagno3no2)!="factor") dt1$flagno3no2<- as.factor(dt1$flagno3no2)
if (class(dt1$flagno2)!="factor") dt1$flagno2<- as.factor(dt1$flagno2)
if (class(dt1$flagnh4)!="factor") dt1$flagnh4<- as.factor(dt1$flagnh4)
if (class(dt1$flagtotnf)!="factor") dt1$flagtotnf<- as.factor(dt1$flagtotnf)
if (class(dt1$flagtotnuf)!="factor") dt1$flagtotnuf<- as.factor(dt1$flagtotnuf)
if (class(dt1$flagtotpf)!="factor") dt1$flagtotpf<- as.factor(dt1$flagtotpf)
if (class(dt1$flagtotpuf)!="factor") dt1$flagtotpuf<- as.factor(dt1$flagtotpuf)
if (class(dt1$flagdrsif)!="factor") dt1$flagdrsif<- as.factor(dt1$flagdrsif)
if (class(dt1$flagbrsif)!="factor") dt1$flagbrsif<- as.factor(dt1$flagbrsif)
if (class(dt1$flagbrsiuf)!="factor") dt1$flagbrsiuf<- as.factor(dt1$flagbrsiuf)
if (class(dt1$flagtpm)!="factor") dt1$flagtpm<- as.factor(dt1$flagtpm)
if (class(dt1$flagtotnuf_sloh)!="factor") dt1$flagtotnuf_sloh<- as.factor(dt1$flagtotnuf_sloh)
if (class(dt1$flagno3no2_sloh)!="factor") dt1$flagno3no2_sloh<- as.factor(dt1$flagno3no2_sloh)
if (class(dt1$flagnh4_sloh)!="factor") dt1$flagnh4_sloh<- as.factor(dt1$flagnh4_sloh)
if (class(dt1$flagkjdl_n_sloh)!="factor") dt1$flagkjdl_n_sloh<- as.factor(dt1$flagkjdl_n_sloh)
if (class(dt1$flagtotpuf_sloh)!="factor") dt1$flagtotpuf_sloh<- as.factor(dt1$flagtotpuf_sloh)
if (class(dt1$flagdrp_sloh)!="factor") dt1$flagdrp_sloh<- as.factor(dt1$flagdrp_sloh)
if (class(dt1$flagdrsif_sloh)!="factor") dt1$flagdrsif_sloh<- as.factor(dt1$flagdrsif_sloh)

# Convert Missing Values to NA for non-dates

# ---- format and export data ----

nuts <- dt1
nuts <- as.data.table(nuts)
nuts[ ,sampledate := parse_date_time(sampledate,"ymd")]

unique(nuts$lakeid)
nuts <- nuts[lakeid == "ME"]

nuts <- nuts[ ,.(year4,daynum,depth,doc,flagdoc,totpuf_sloh,flagtotpuf_sloh,drp_sloh,flagdrp_sloh)]

# Total Phosphorus

tp <- nuts[depth <= 12 ,.(year4,daynum,depth,totpuf_sloh, flagtotpuf_sloh)]

# remove crazy values
unique(tp$flagtotpuf_sloh) # remove J and Q flags, R and S seem OK, just low value flags
tp <- tp[flagtotpuf_sloh != "Q" & flagtotpuf_sloh != "J"]
tp <- tp[!is.na(totpuf_sloh)]
# nothing extreme

tp[ , .N, by = .(year4,depth)]
# I dunno what's best, seems problematic to interpolate by depth to the lower depths when they're missing
# But also don't want to just use the surface measurements
# Take average of all available epilimnion measurements on each day
# Then take average of all days
# These are pretty evenly spaced, so again I think it's OK to not interpolate

tp <- tp[ , .(Mean.TP.mg.L = mean(totpuf_sloh, na.rm = TRUE)), by = .(year4, daynum)] # first average across available epi measurements on each day
tp <- tp[ , .(Mean.TP.mg.L = mean(Mean.TP.mg.L, na.rm = TRUE)), by = .(year4)] # next average across years

tp <- tp[!is.na(Mean.TP.mg.L)]

# Soluble Reactive Phosphorus

srp <- nuts[depth <= 12 ,.(year4,daynum,depth,drp_sloh, flagdrp_sloh)]

# remove crazy values
unique(srp$flagdrp_sloh)
srp <- srp[flagdrp_sloh != "Q" & flagdrp_sloh != "J"]
srp <- srp[!is.na(drp_sloh)]
# nothing extreme

srp <- srp[ , .(Mean.SRP.mg.L = mean(drp_sloh, na.rm = TRUE)), by = .(year4, daynum)] # first average across available epi measurements on each day
srp <- srp[ , .(Mean.SRP.mg.L = mean(Mean.SRP.mg.L, na.rm = TRUE)), by = .(year4)] # next average across years

srp <- srp[!is.na(Mean.SRP.mg.L)]

# DOC

doc <- nuts[depth <= 12 ,.(year4,daynum,depth,doc, flagdoc)]

# remove crazy values
unique(doc$flagdoc)
doc <- doc[flagdoc != "L" & flagdoc != "K" & flagdoc != "J" & flagdoc != "H" & flagdoc != "G" & flagdoc != "GK" & flagdoc != "GEK" & flagdoc != "FGK"]
doc <- doc[!is.na(doc)]
doc <- doc[doc >= 0]

doc <- doc[ , .(Mean.DOC.mg.L = mean(doc, na.rm = TRUE)), by = .(year4, daynum)] # first average across available epi measurements on each day
doc <- doc[ , .(Mean.DOC.mg.L = mean(Mean.DOC.mg.L, na.rm = TRUE)), by = .(year4)] # next average across years

doc <- doc[!is.na(Mean.DOC.mg.L)]

# ---- Export data ----

fwrite(x = tp, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/total_phosphorus.csv")
fwrite(x = srp, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/soluble_reactive_phosphorus.csv")
fwrite(x = doc, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/dissolved_organic_carbon.csv")
