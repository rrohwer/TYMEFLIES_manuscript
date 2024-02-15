# RRR
# 6/9/21
# I'm going to re-name these unintelligible files from jamo
# make a reference of filename, spID, and my name

jamo <- read.delim(file = "data/2021-06-07_jamo_fetching_spids/2021-06-09_spID_status.txt", header = F, sep = " ", colClasses = "character")

my.ref <- readRDS(file = "data/2021-05-24_metaG_stats/all_img_assemblies.rds")

created.file <- "data/2021-06-07_jamo_fetching_spids/2021-06-09_key_with_jamofile_spid_tymename.txt"

# ----

# jamo is the output of "jamo info spid 1263380 >> spID_status.txt" repeated for each spid

head(jamo)
jamo <- jamo[ ,-5]
colnames(jamo) <- c("ITS.SP.ID","jamo.filename","jamo.status.6.9.21","jamo.file.ID")
unique(jamo$jamo.status.6.9.21) # all restored from tape
head(jamo)
str(jamo)

head(my.ref)
str(my.ref)

combo <- merge(x = jamo, y = my.ref, by = "ITS.SP.ID", all.x = T, all.y = F)
head(combo)
colnames(combo)

combo <- combo[ ,c(1,2,4,11)]
head(combo)

# ----

created.file
# write.table(x = combo, file = created.file, sep = "\t", quote = F, row.names = F)
