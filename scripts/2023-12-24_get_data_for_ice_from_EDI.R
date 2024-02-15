# RRR

library(data.table)
library(lubridate)

# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.33.37

# ---- EDI download ----

# Package ID: knb-lter-ntl.33.37 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Ice Duration - Madison Lakes Area 1853 - current.
# Data set creator:  John Magnuson - University of Wisconsin-Madison 
# Data set creator:  Stephen Carpenter - University of Wisconsin-Madison 
# Data set creator:  Emily Stanley - University of Wisconsin-Madison 
# Metadata Provider:  NTL Information Manager - University of Wisconsin-Madison 
# Contact:    -  NTL LTER  - ntl.infomgr@gmail.com
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/33/37/933b6eb9b31cc1e41c6a02ee40a91877" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lakeid",     
                 "season",     
                 "iceon",     
                 "ice_on",     
                 "iceoff",     
                 "ice_off",     
                 "ice_duration",     
                 "year4"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$lakeid)!="factor") dt1$lakeid<- as.factor(dt1$lakeid)
if (class(dt1$season)!="factor") dt1$season<- as.factor(dt1$season)
if (class(dt1$iceon)=="factor") dt1$iceon <-as.numeric(levels(dt1$iceon))[as.integer(dt1$iceon) ]               
if (class(dt1$iceon)=="character") dt1$iceon <-as.numeric(dt1$iceon)                                   
# attempting to convert dt1$ice_on dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1ice_on<-as.Date(dt1$ice_on,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1ice_on) == length(tmp1ice_on[!is.na(tmp1ice_on)])){dt1$ice_on <- tmp1ice_on } else {print("Date conversion failed for dt1$ice_on. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1ice_on) 
if (class(dt1$iceoff)=="factor") dt1$iceoff <-as.numeric(levels(dt1$iceoff))[as.integer(dt1$iceoff) ]               
if (class(dt1$iceoff)=="character") dt1$iceoff <-as.numeric(dt1$iceoff)                                   
# attempting to convert dt1$ice_off dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1ice_off<-as.Date(dt1$ice_off,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1ice_off) == length(tmp1ice_off[!is.na(tmp1ice_off)])){dt1$ice_off <- tmp1ice_off } else {print("Date conversion failed for dt1$ice_off. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1ice_off) 
if (class(dt1$ice_duration)=="factor") dt1$ice_duration <-as.numeric(levels(dt1$ice_duration))[as.integer(dt1$ice_duration) ]               
if (class(dt1$ice_duration)=="character") dt1$ice_duration <-as.numeric(dt1$ice_duration)
if (class(dt1$year4)=="factor") dt1$year4 <-as.numeric(levels(dt1$year4))[as.integer(dt1$year4) ]               
if (class(dt1$year4)=="character") dt1$year4 <-as.numeric(dt1$year4)

# Convert Missing Values to NA for non-dates

dt1$iceon <- ifelse((trimws(as.character(dt1$iceon))==trimws("NA")),NA,dt1$iceon)               
suppressWarnings(dt1$iceon <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$iceon))==as.character(as.numeric("NA"))),NA,dt1$iceon))
dt1$iceoff <- ifelse((trimws(as.character(dt1$iceoff))==trimws("NA")),NA,dt1$iceoff)               
suppressWarnings(dt1$iceoff <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$iceoff))==as.character(as.numeric("NA"))),NA,dt1$iceoff))
dt1$ice_duration <- ifelse((trimws(as.character(dt1$ice_duration))==trimws("NA")),NA,dt1$ice_duration)               
suppressWarnings(dt1$ice_duration <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$ice_duration))==as.character(as.numeric("NA"))),NA,dt1$ice_duration))
dt1$year4 <- ifelse((trimws(as.character(dt1$year4))==trimws("NA")),NA,dt1$year4)               
suppressWarnings(dt1$year4 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$year4))==as.character(as.numeric("NA"))),NA,dt1$year4))


# ---- save data for plot ----              

ice <- as.data.table(dt1)
ice <- ice[lakeid == "ME"]
ice[ ,summer.year := year4 + 1]
ice <- ice[ ,.(summer.year, ice_duration)]

fwrite(x = ice, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/ice_durations.csv")
