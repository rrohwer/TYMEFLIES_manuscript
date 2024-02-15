# RRR
# Will run this from dtn (data transfer node)

# path.file.info <- "finding-data/2021-12-23_filtered_reads/uniqueSPIDS.txt"
path.file.info <- "finding-data/2021-12-23_filtered_reads/dupSPIDs_made_unique.txt"
path.new.folder <- "/global/cfs/cdirs/pkmeco/tymeflies/filtered_reads/"


file.info <- read.delim(file = path.file.info, header = T, sep = " ")
head(file.info)

my.command <- "cp -i"
old.location <- file.info$FILENAME
new.location <- path.new.folder

new.script <- paste(my.command,old.location,new.location)

# write.table(x = new.script, row.names = F, quote = F, col.names = F, 
#             file = "finding-data/2021-12-23_filtered_reads/2021-12-23_3_move_files.sh")

write.table(x = new.script, row.names = F, quote = F, col.names = F, 
            file = "finding-data/2021-12-23_filtered_reads/2022-01-20_6_move_files-dupSPIDS.sh")
