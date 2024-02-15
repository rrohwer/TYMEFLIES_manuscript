# RRR
# this should not be slow, don't bother with parallel. Just read in the combined file, the genome file for coverage, and the tax file
# and then make the plots one after the other. 
# maybe even just run on my laptop?

library(data.table)
library(lubridate)

file.selection <- "data/2023-10-26_selection_per_sample_summaries/all_genome_selection_summaries.tsv.gz"
file.genome <- "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz"
file.tax <- "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds"

breadth.cutoff <- .5
coverage.cutoff <- 10

# ---- 

selection <- fread(file = file.selection)
selection <- selection[ ,.(genome, sample, Npos, Nneg)]

genome <- fread(file = file.genome)
genome[ ,date := parse_date_time(x = date, orders = "ymd")]
genome[ ,above.breadth.cutoff := (breadth / breadth_expected) >= breadth.cutoff, ]
genome <- genome[above.breadth.cutoff == TRUE & coverage_median > coverage.cutoff]

tax <- readRDS(file = file.tax)
tax <- as.data.table(tax)
tax <- tax[winner == TRUE]

for (g in tax$bin.full.name){
  my.genome <- genome[genome == g]
  
  if (nrow(my.genome) > 0){
    my.genome <- merge(x = my.genome, y = selection, by = c("genome","sample"), all.x = TRUE, all.y = FALSE)
    my.genome[is.na(Npos) ,`:=`(Npos = 0)]
    my.genome[is.na(Nneg) ,`:=`(Nneg = 0)]
    
    my.tax <- tax[bin.full.name == g]
    
    by.year <- my.genome[ ,.(.N, 
                             Npos.tot = sum(Npos), 
                             Nneg.tot = sum(Nneg),
                             Npos.ave = mean(Npos),
                             Npos.sd = sd(Npos), 
                             Nneg.tot = sum(Nneg), 
                             Nneg.tot = sum(Nneg),
                             Nneg.ave = mean(Nneg),
                             Nneg.sd = sd(Nneg)), by = year]
    
    by.season <- my.genome[ ,.(.N, 
                             Npos.tot = sum(Npos), 
                             Nneg.tot = sum(Nneg),
                             Npos.ave = mean(Npos),
                             Npos.sd = sd(Npos), 
                             Nneg.tot = sum(Nneg), 
                             Nneg.tot = sum(Nneg),
                             Nneg.ave = mean(Nneg),
                             Nneg.sd = sd(Nneg)), by = season]
    
    # make plot of all samples individually
    x.axis
    plot(x = my.genome[ ,date], y = my.genome[ ,Npos], ylim = c(-max(my.genome[ ,Nneg]), max(my.genome[ ,Npos])), type = "n")
    
  }
  
}