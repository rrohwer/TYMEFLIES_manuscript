# RRR

library(readxl)

folder.on.cfs <- "/global/cfs/cdirs/pkmeco/tymeflies/assembled_contigs"

path.bigtyme <- "data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx"

# Choose:

file.spid.info <- "finding-data/2022-03-08_assembled_contigs/uniqueSPIDS.txt"
created.filename.key <- "finding-data/2022-03-08_assembled_contigs/filename_key_unique.rds"
created.shell.script <- "finding-data/2022-03-08_assembled_contigs/4a-get_contigs-copy_and_rename-uniqueSPIDs.sh"

# file.spid.info <- "finding-data/2022-03-08_assembled_contigs/dupSPIDS_chosen.txt"
# created.filename.key <- "finding-data/2022-03-08_assembled_contigs/filename_key_dupSPIDs.rds"
# created.shell.script <- "finding-data/2022-03-08_assembled_contigs/4b-get_contigs-copy_and_rename-dupSPIDs.sh"


# get existing file names in jamo ----

jamo.filenames <- read.delim(file = file.spid.info, header = T, sep = " ")
head(jamo.filenames)
tail(jamo.filenames)

# create new unique names for the files ----

bigtyme <- read_excel(path = path.bigtyme, na = "NA")
bigtyme[1:5,1:5]
colnames(bigtyme)

big <- merge(x = bigtyme, y = jamo.filenames, by.x = "ITS_SPID", by.y = "SPID")
head(big)
big$FILENAME
index.missing <-which(big$FILENAME == "" | is.na(big$FILENAME) )
if(length(index.missing) > 0){
  big$ITS_SPID[index.missing]
  big <- big[-index.missing, ]
}

big$new.names <- paste0(big$TYMEFLIES.name,"_",big$taxon_oid,"_","assembled_contigs.fasta" )

# save a filename key ----

key.filenames <- big[ ,c(1,39,40)]
head(key.filenames)

# create a shell script to copy and rename the files

my.command <- "cp -i"
old.location <- big$FILENAME
new.location <- paste0(folder.on.cfs,"/",big$new.names)

new.script <- paste(my.command,old.location,new.location)
head(new.script)

# export files ----

cat(created.filename.key)
saveRDS(object = key.filenames, file = created.filename.key)

cat(created.shell.script)
write.table(x = new.script, row.names = F, quote = F, col.names = F, file = created.shell.script)

# end ----