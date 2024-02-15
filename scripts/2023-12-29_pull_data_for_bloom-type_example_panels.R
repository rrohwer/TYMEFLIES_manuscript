# RRR 
# abundance panel of the less diverse blooms example plot
# fig 2c,2d
# example genomes: 
# less diverse bloom/anicorrelated: ME2011-09-04_3300044729_group3_bin142
# more diverse bloom/correlated: ME2012-08-31_3300044613_group4_bin150

library(data.table)
library(lubridate)

genome.info <- fread(file = "data/2023-12-06_abundance_seasonality_analysis/filtered_abundance_and_nucl_div.tsv.gz", colClasses = c("date" = "character"))

genome.info[ ,date := parse_date_time(date, "ymd")]
genome.info[ ,`:=`(year = year(date), yday = yday(date))]
genome.info[ ,date := as.character(date)]
genome.info <- genome.info[order(date)]

clonal.example <- genome.info[genome == "ME2011-09-04_3300044729_group3_bin142", ]
diverse.example <- genome.info[genome == "ME2012-08-31_3300044613_group4_bin150"]

# calculate moving average to show a smoothed line in figure
window <- 30
points.freq <- 1

# clonal bloom
moving.ave <- data.table("floor" = seq(from = 0, to = 365, by = points.freq))
moving.ave[ ,ceiling := floor + window]
moving.ave[ ,center := floor + (window / 2)]
for(r in 1:nrow(moving.ave)){
  my.vals <- clonal.example[yday >= moving.ave[r,floor] & yday <= moving.ave[r,ceiling], nucl_diversity]
  moving.ave[r, ave.div := mean(my.vals, na.rm = T)]
  moving.ave[r, sd.div := sd(my.vals, na.rm = T)]
  moving.ave[r, num.obs.div := length(my.vals[!is.na(my.vals)])]
  my.vals <- clonal.example[yday >= moving.ave[r,floor] & yday <= moving.ave[r,ceiling], abund.perc]
  moving.ave[r, ave.abund := mean(my.vals)]
  moving.ave[r, sd.abund := sd(my.vals)]
  moving.ave[r, num.obs.abund := length(my.vals)]
}
moving.ave.clonal <- moving.ave

# diverse bloom
moving.ave <- data.table("floor" = seq(from = 0, to = 365, by = points.freq))
moving.ave[ ,ceiling := floor + window]
moving.ave[ ,center := floor + (window / 2)]
for(r in 1:nrow(moving.ave)){
  my.vals <- diverse.example[yday >= moving.ave[r,floor] & yday <= moving.ave[r,ceiling], nucl_diversity]
  moving.ave[r, ave.div := mean(my.vals, na.rm = T)]
  moving.ave[r, sd.div := sd(my.vals, na.rm = T)]
  moving.ave[r, num.obs.div := length(my.vals[!is.na(my.vals)])]
  my.vals <- diverse.example[yday >= moving.ave[r,floor] & yday <= moving.ave[r,ceiling], abund.perc]
  moving.ave[r, ave.abund := mean(my.vals)]
  moving.ave[r, sd.abund := sd(my.vals)]
  moving.ave[r, num.obs.abund := length(my.vals)]
}
moving.ave.diverse <- moving.ave

# ---- save data ----

fwrite(x = clonal.example, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/less_diverse_bloom-ME2011-09-04_3300044729_group3_bin142.csv")
fwrite(x = diverse.example, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/more_diverse_bloom-ME2012-08-31_3300044613_group4_bin150.csv")

fwrite(x = moving.ave.clonal, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/less_diverse_bloom_1mo_moving_average-ME2011-09-04_3300044729_group3_bin142.csv")
fwrite(x = moving.ave.diverse, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_2_data/more_diverse_bloom_1mo_moving_average-ME2012-08-31_3300044613_group4_bin150.csv")

