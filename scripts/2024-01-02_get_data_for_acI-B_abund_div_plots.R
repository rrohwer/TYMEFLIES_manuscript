# RRR
# data for the acI-B figure: ME2011-09-21_3300043464_group3_bin69


library(data.table)

genome <- fread(file = "data/2023-12-08_combined_genome_info/combined_genome_info_by_sample.tsv.gz", colClasses = c("date" = "character"))

genome <- genome[genome == "ME2011-09-21_3300043464_group3_bin69"]

abund.div <- genome[ ,.(date, year, adj.abund.perc, nucl_diversity)]

color.key <- data.table("year" = 2000:2019, 
                        "color.year.step" = c(rep("#b89996", 11), "purple","red","tan2","dodgerblue4", rep("#7592b2", 5)))

abund.div <- merge(x = abund.div, y = color.key, by = "year")


fwrite(x = abund.div, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/abund_and_nucl_div-ME2011-09-21_3300043464_group3_bin69.csv")
