# RRR
# pull file from the server:
# multivariate_SNV_analysis-med_cov_10/NMDS_objects_all_SNVs/ME2011-09-21_3300043464_group3_bin69_all_SNV_euclidean_distance_NMDS_object.rds

library(data.table)
library(lubridate)

sample.key <- fread(file = "data/2023-12-08_combined_genome_info/sample_key.tsv", colClasses = c("date" = "character"))

my.nmds <- readRDS(file = "data/2023-12-29_genome_examples_data_for_paper_figures/ME2011-09-21_3300043464_group3_bin69_all_SNV_euclidean_distance_NMDS_object.rds")

# ---- set up for plot ----

sample.key[ ,date := parse_date_time(date,"ymd")]

my.nmds <- data.table("sample" = rownames(my.nmds$points),"x" = my.nmds$points[ ,1], "y" = my.nmds$points[ ,2])
my.nmds <- merge(x = my.nmds, y = sample.key, by = "sample")
my.nmds <- my.nmds[order(date)]

my.nmds <- my.nmds[ ,.(date,Year,x,y)]

color.key <- data.table("Year" = 2000:2019, 
                        "color.year.step" = c(rep("#b89996", 11), "purple","red2","tan2","dodgerblue4", rep("#7592b2", 5)))
my.nmds <- merge(x = my.nmds, y = color.key, by = "Year")               

my.nmds[ ,date := as.character(date)]

fwrite(x = my.nmds, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/nmds_coords-ME2011-09-21_3300043464_group3_bin69.csv")



