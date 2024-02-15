# RRR
# data for the acI-B figure: ME2011-09-21_3300043464_group3_bin69

library(data.table)

selection <- fread(file = "data/2023-10-26_selection_per_sample_summaries/all_genome_selection_summaries.tsv.gz", colClasses = c("date" = "character"))

selection <- selection[genome == "ME2011-09-21_3300043464_group3_bin69"]

color.key <- data.table("year" = 2000:2019, 
                        "color.year.step" = c(rep("#b89996", 11), "purple","red2","tan2","dodgerblue4", rep("#7592b2", 5)))

selection <- merge(x = selection, y = color.key, by = "year")

fwrite(x = selection, file = "figures/2023-12-24_paper_figures/Rohwer_Figure_5_data/amt_selection-ME2011-09-21_3300043464_group3_bin69.csv")


