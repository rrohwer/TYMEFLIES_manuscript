# RRR

library(data.table)
library(lubridate)
library(limony)

data("seasons") # from limony

# get the start/stop season dates for falls and ice-ons that span the new year ----

seasons # these are season start dates, need to add more ice details for plotting
seasons <- as.data.table(seasons)

seasons[year(`Ice-on`) == Year, fall.ice.on := `Ice-on`]
seasons[is.na(fall.ice.on), fall.ice.on := parse_date_time(paste(Year, 12, 31), "ymd", tz = "Etc/GMT-5")] 
seasons[ ,fall.ice.off := parse_date_time(paste(Year, 12, 31), "ymd", tz = "Etc/GMT-5")]

seasons[ , prev.ice.on := c(seasons$`Ice-on`[1],seasons$`Ice-on`[-nrow(seasons)])] 
seasons[year(prev.ice.on) == Year, spring.ice.on := prev.ice.on]
seasons[year(prev.ice.on) == Year - 1, spring.ice.on := parse_date_time(paste(Year, 1, 1), "ymd", tz = "Etc/GMT-5")]

seasons[year(`Ice-on`) == Year, fall.fall.end := `Ice-on`]
seasons[year(`Ice-on`) == Year + 1, fall.fall.end := parse_date_time(paste(Year, 12, 31), "ymd", tz = "Etc/GMT-5")]

years.where.fall.crosses.new.year <- which(yday(seasons$spring.ice.on) > 1)
seasons[years.where.fall.crosses.new.year, spring.fall.start := parse_date_time(paste(Year, 1, 1), "ymd", tz = "Etc/GMT-5")]

seasons <- seasons[Year > 1999]
seasons <- seasons[ ,-c("prev.ice.on")]

# put in sequential order ----
cbind(1:ncol(seasons), colnames(seasons))
seasons <- seasons[ ,c(1,12,13,10,2,3,4,5,6,11,8,9)]

# save this table and recover from headache ----

saveRDS(object = seasons, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_1_data/season_dates.rds")

