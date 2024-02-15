# RRR

# have the GFF annotation files
# now just need to add that to the gene info files 
# and hope they match!

# first we have to do step 5- concatenate all the gff files into one!
# $ cat ../../yggshare/current_projects/TYMEFLIES/tymeflies/gene_annotations_IMG/*.gff > all_IMG_annotations.gff
# $ ls -lh
# -rw-rw-r-- 1 rrohwer rrohwer 178G Apr 12 20:41 all_IMG_annotations.gff
# it's a pretty big file!

# NOTE:
# instrain gene start is zero-indexed
# instrain gene end is zero-indexed
# gff start has a 1-base offset
# gff end has a 1-base offset
# I think this means I subtract 1 from every start/stop position in the gff file


library(data.table)

userprefs <- commandArgs(trailingOnly = TRUE)
key <- userprefs[1]
genes <- userprefs[2]
gff <- userprefs[3]
output <- userprefs[4]
threads <- userprefs[5]

threads <- as.numeric(threads)

setDTthreads(threads = threads)

# cat("you forgot to comment out the file paths!")
# key <- fread(file = "data/2023-04-12_get_contig_name_key_testing/small_bin_scaffold_key.csv")
# genes <- fread(file = "data/2023-04-12_get_contig_name_key_testing/small_gene_info_combined.tsv")
# gff <- fread(file = "data/2023-04-12_get_contig_name_key_testing/small_gene_annotation.gff")

cat("reading in the scaffold key file.\n")
key <- fread(file = key, nThread = threads, verbose = T)

cat("reading in the instrain gene info file. \n")
genes <- fread(file = genes, nThread = threads, verbose = T)

cat("reading in the huuuge gene annotation gff file.\n")
gff <- fread(file = gff, nThread = threads, verbose = T)

cat("minor formatting.\n")
colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
key[ ,Assembly.Scaffold.Name := NULL]

cat("shifting the start/stop in gff to be 0-based.\n")
gff[ , `:=`(start = start - 1, end = end - 1)]

cat("merging the gff file with the scaffold key file.\n")
gff <- merge(x = gff, y = key, by.x = "seqid", by.y = "IMG.Scaffold.Name", all.x = FALSE, all.y = FALSE)

cat("merging the gff file with the instrain gene info file.\n")
genes <- merge(x = genes, y = gff, by.x = c("scaffold","start","end"), by.y = c("Robin.Scaffold.Name","start","end"), all.x = TRUE, all.y = FALSE)

cat("saving the new gene info file.\n")
fwrite(x = genes, file = output, sep = "\t", nThread = threads)



