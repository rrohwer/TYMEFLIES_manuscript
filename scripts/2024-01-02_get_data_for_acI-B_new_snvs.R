# RRR
# data for the acI-B figure: ME2011-09-21_3300043464_group3_bin69
# pulled data off of server: new_SNV_analysis-med_cov_10/new_SNV_stats/ME2011-09-21_3300043464_group3_bin69_new_and_total_SNVs.tsv.gz

library(data.table)

snv <- fread(file = "data/2023-12-29_genome_examples_data_for_paper_figures/ME2011-09-21_3300043464_group3_bin69_new_and_total_SNVs.tsv.gz", colClasses = c("Date" = "character"))

snv <- snv[!is.na(New) ,.(Date, year, New)]

color.key <- data.table("year" = 2000:2019, 
                        "color.year.step" = c(rep("#b89996", 11), "purple","red2","tan2","dodgerblue4", rep("#7592b2", 5)))

snv <- merge(x = snv, y = color.key, by = "year")

fwrite(x = snv, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/new_snvs_ME2011-09-21_3300043464_group3_bin69.csv")


# # check
# snv[ ,Date := parse_date_time(Date, "ymd")]
# 
# plot(snv$Date[-c(1:10)], snv$New[-c(1:10)], type = "l") # exclude at least first 9
