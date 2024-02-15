if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# It gave an error message that it couldn't install the dependency "genefilter" package
# solution- install the intel-version of R. Bioconductor can't do M1 mac yet. sigh

browseVignettes("DESeq2") # damnit not installed

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefilter")

library(DESeq2)

