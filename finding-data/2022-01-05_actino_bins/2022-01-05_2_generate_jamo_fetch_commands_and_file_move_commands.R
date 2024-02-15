# RRR
# so now I parse the jamo info output into jamo fetch commands
# for actino IMG bins containing 16S
# notice my grep command added somefiles, I think grep 33000123_2 also returns 33000123_22 for example

jamo.out <- read.delim(file = "finding-data/2022-01-05_actino_bins/2022-01-05_1_jamo_info.txt", header = F, sep = " ")

bins <- read.csv(file = "data/2022-01-05_IMG_IDs_for_Actino_bins/IMG_bins_Actinobacteria_with16S.csv")

created.fetch.script <- "finding-data/2022-01-05_actino_bins/2022-01-05_2_fetch_bins_jamo.sh"
created.move.script <- "finding-data/2022-01-05_actino_bins/2022-01-05_2_copy_bins_into_directory.sh"

head(jamo.out)
head(bins)

jamo.out <- jamo.out[ ,-5]
colnames(jamo.out) <- c("ITS_SPID","FILENAME","STATUS","FILEID")
head(jamo.out)

# notice my grep command added some files, I think grep 33000123_2 also returns 33000123_22 for example
index <- duplicated(jamo.out$FILENAME) # Oh but also the exact file is duplicated like it got returned twice
jamo.out <- jamo.out[!index, ]

# now also remove the bins where different numbers matched the grep
jamo.out$Bin.ID <- sub(pattern = "^.*/", replacement = "", x = jamo.out$FILENAME)
jamo.out$Bin.ID <- sub(pattern = "\\.tar\\.gz$", replacement = "", x = jamo.out$Bin.ID)

# wait there's also missing bins. checked and they are not in jamo
all(bins$Bin.ID %in% jamo.out$Bin.ID) # why not???? But this is real, not a bug in this script. looked individually in jamo.
index <- which(bins$Bin.ID %in% jamo.out$Bin.ID)
missing.bins <- bins[-index, ]
write.table(x = missing.bins, file = "~/Desktop/bins_missing_from_jamo.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

# get the bins I could find
new.script.fetch <- paste("jamo fetch all id", jamo.out$FILEID)

new.script.copy <- paste("cp", jamo.out$FILENAME, "/global/cscratch1/sd/rohwer/img_actino_bins_with_16S/")

# ---- export ----

write.table(x = new.script.fetch, file = created.fetch.script, quote = F, row.names = F, col.names = F)
write.table(x = new.script.copy, file = created.move.script, quote = F, row.names = F, col.names = F)

# end