# RRR
# Will run this from dtn (data transfer node)

path.file.info <- "finding-data/2021-12-23_filtered_reads/uniqueSPIDS.txt"
path.new.folder <- "/global/cfs/cdirs/pkmeco/tymeflies/filtered_reads/"


file.info <- read.delim(file = path.file.info, header = T, sep = " ")
head(file.info)

for (f in file.info$FILENAME){
  old.location <- f
  new.location <- path.new.folder
  
  args.str <- paste("-i",old.location,new.location)
}

system2(command = "cp", args = args.str, stdout = "termout_2021-12-23.txt")
