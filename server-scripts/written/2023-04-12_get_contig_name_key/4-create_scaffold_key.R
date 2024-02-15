# RRR
# make scaffold key
# will need to run on server b/c of large size of the IMG key that has alll scaffolds in it
# will create a key with just the bin scaffolds, since those are the only ones I renamed.

# ---- set up ----

library(data.table)

userprefs <- commandArgs(trailingOnly = TRUE)
file.bigtyme <- userprefs[1]
file.bin.scaffolds <- userprefs[2]
file.img.scaffolds <- userprefs[3]
output.file <- userprefs[4]
threads <- userprefs[5]

# # for testing
# cat("\nyou forgot to re-comment file paths!\n")
# file.bigtyme <- "data/2022-06-27_binning_groups/2022-06-27_bigtyme.rds"
# file.bin.scaffolds <- "data/2023-04-12_get_contig_name_key_testing/small_bin_scaffolds.txt"
# file.img.scaffolds <- "data/2023-04-12_get_contig_name_key_testing/small_IMG_scaffolds.txt"
# output.file <- "data/2023-04-12_get_contig_name_key_testing/test_output.txt"

# ---- read in files ----

threads <- as.numeric(threads)
bigtyme <- readRDS(file = file.bigtyme)
bin.scaffolds <- fread(file = file.bin.scaffolds, header = F, nThread = threads)
img.scaffolds <- fread(file = file.img.scaffolds, header = F, nThread = threads)

# ---- parse out each file ----

bigtyme <- as.data.table(bigtyme)
colnames(bigtyme)[colnames(bigtyme) == "GOLD Analysis Project ID"] <- "GOLD.AP.ID"
bigtyme <- bigtyme[ ,.(Sample.Name = sample.names, GOLD.AP.ID)]

bin.scaffolds$V1 <- sub("^>", "", bin.scaffolds$V1)
bin.scaffolds <- bin.scaffolds[ ,.(Robin.Scaffold.Name = V1)]
bin.scaffolds[ ,Sample.Name := sub(pattern = "_group.*", replacement = "", x = Robin.Scaffold.Name)]
bin.scaffolds[ ,Assembly.Scaffold.Name := sub(pattern = "^.*\\.fna_", replacement = "", x = Robin.Scaffold.Name)]

bin.scaffolds <- merge(x = bin.scaffolds, y = bigtyme, by = "Sample.Name", all.x = TRUE, all.y = FALSE)

img.scaffolds <- img.scaffolds[ ,.(Assembly.Scaffold.Name = V1, IMG.Scaffold.Name = V2)]
img.scaffolds[ ,GOLD.AP.ID := sub(pattern = "_.*$", replacement = "", x = IMG.Scaffold.Name)]

img.scaffolds <- merge(x = img.scaffolds, y = bin.scaffolds, by = c("Assembly.Scaffold.Name", "GOLD.AP.ID"), all.x = FALSE, all.y = TRUE) 

img.scaffolds <- img.scaffolds[ ,.(Assembly.Scaffold.Name,IMG.Scaffold.Name,Robin.Scaffold.Name)]

# ---- save resulting key ----

fwrite(x = img.scaffolds, file = output.file, nThread = threads) # makes csv



