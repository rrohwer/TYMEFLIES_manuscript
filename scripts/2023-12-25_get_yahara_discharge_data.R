# RRR
# download the latest yahara discharge data from USGS
# get it formatted for the 2012 plot

library(dataRetrieval) # USGS data download package
library(data.table)
library(lubridate)

# ---- import ----

# Yahara at Windsor is the yaraha input location that goes back the longest
# site ID: 05427718
# discharge parameter ID: 00060
# P loading parameter ID: 91050
# N loading (TN) param ID: 00601
# find these codes here: https://help.waterdata.usgs.gov/codes-and-parameters/parameters
# also need that ^ to figure out the units!
# package help here: https://doi-usgs.github.io/dataRetrieval/
# (not very helpful)
# and how to cite:
# https://help.waterdata.usgs.gov/faq/miscellaneous/how-to-cite-usgs-water-data-for-the-nation-waterdata.usgs.gov-in-a-publication

discharge <- readNWISdv(siteNumbers = "05427718", parameterCd = "00060")
discharge <- as.data.table(discharge)
colnames(discharge)[4] <- "discharge.ft3.s" # cubic feet per second

# ---- format and export ----

discharge[ ,year := year(Date)]
discharge[ ,yday := yday(Date)]

discharge <- discharge[ ,c(1,2,3,6,7,4)]

discharge[ ,Date := as.character(Date)]

fwrite(x = discharge, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_4_data/yahara_discharge.csv")
