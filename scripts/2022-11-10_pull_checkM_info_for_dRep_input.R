
input.data <- "data/2022-09-23_bin_stats/all_bins.rds"
output.data.folder <- "data/2022-11-10_dRep_files/"


# ---- get input checkM results ----

bins <- readRDS(input.data)
colnames(bins)

bins <- bins[ ,1:3]
colnames(bins) <- c("genome", "completeness", "contamination")

bins$genome <- paste0(bins$genome,".fna")

write.csv(x = bins, file = file.path(output.data.folder,"bin_info_for_dRep.csv"), quote = F, row.names = F)


# ---- get input bin locations file ----

folder.names <- sub(pattern = "_group.*_bin.*$", replacement = "", x = bins$genome)

bin.locs <- paste("metabat2_bins_renamed",folder.names, bins$genome, sep = "/")

write.table(x = bin.locs, file = file.path(output.data.folder,"bin_locations_for_dRep.txt"), quote = F, row.names = F, col.names = F)


# ---- get test files ----


summary(bins)
# one contamination is >100, may cause errors with dRep

# get subset to test dRep
index.sample1 <- grep(pattern = "ME2000-03-15", x = bins$genome, value = F)
index.sample2 <- grep(pattern = "ME2000-03-30", x = bins$genome, value = F)

test.bins <- bins[c(index.sample1,index.sample2), ]

test.locs <- bin.locs[c(index.sample1,index.sample2)]

write.csv(x = test.bins, file = file.path(output.data.folder,"test_bin_info_for_dRep.csv"), quote = F, row.names = F)
write.table(x = test.locs, file = file.path(output.data.folder,"test_bin_locations_for_dRep.txt"), quote = F, row.names = F, col.names = F)
