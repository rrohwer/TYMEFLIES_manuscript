# RRR

# this package is A. in R, B. actively maintained and cited, and C. has over 1000 citations total

# BUT, it seems euk oriented

# https://yulab-smu.top/biomedical-knowledge-mining-book/

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
  BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("GOSemSim")

library(GOSemSim) # GOSemSim v2.24.0


# OK I don't actually think this does what I want.
# it is set up to compare individual GO terms, not to summarize a whole table of GO terms

# fuck GO terms. seems like nobody uses them or IPR accessions in metaG analyses

# move on to KEGG