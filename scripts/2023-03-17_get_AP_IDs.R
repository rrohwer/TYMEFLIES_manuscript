# RRR
# get AP IDs for Matthew to download for me.

x <- readRDS("data/2022-06-27_binning_groups/2022-06-27_bigtyme.rds")
colnames(x)
unique(x$NOTES)
apid <- x$`ITS AP ID`
write.table(x = apid, file = "~/Desktop/TYMEFLIES_AP_IDs.txt", quote = F, row.names = F, col.names = F)
