# RRR

# parse the genome names to prep for kegg decoder before combining all the kegg annotation files
# this script will be for one genome at a time, and I'll run them in parallel

library(data.table)

userinput <- commandArgs(trailingOnly = TRUE)

input <- userinput[1]
output <- userinput[2]

# input <- "data/2023-09-13_KEGG_annotations/ME2011-09-21_3300043464_group3_bin69.kofamscan.tsv"
# output <- "data/2023-09-13_KEGG_annotations/ME2011-09-21_3300043464_group3_bin69.kofamscan-no_underscores.tsv"

ko <- fread(file = input, sep2 = "\t", header = T, nThread = 1)

ko <- ko[ ,-1]
ko[ ,`gene name` := sub(pattern = "_", replacement = ".", x = `gene name`)]
ko[ ,`gene name` := sub(pattern = "_", replacement = ".", x = `gene name`)]
ko[ ,`gene name` := sub(pattern = "_", replacement = ".", x = `gene name`)]

fwrite(x = ko, file = output, sep = "\t")
