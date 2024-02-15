# RRR
# similar to coverM genome-definition file , inStrain needs same file with cols reversed
# make for both coverM and inStrain here:

# coverM wants:
# --genome-definition		FILE containing list of genome_name<tab>contig lines to define the genome of each contig.

# inStrain wants:
# scaffold-to-bin file  text file with two columns separated by tabs, where the first column is the name of a scaffold and the second column is the name of the bin / genome the scaffold belongs to. 

# ---- get list of scaffold names ----

# I pulled out the list of contigs this way
# cat drep_winners_96_concat.fna | grep ">" > drep_winners_96_contig_list.txt
# read that in and parse it into the input file I need

# ---- make input files ----

contigs <- read.delim(file = "data/2023-03-13_inStrain_on_drep96/input/drep_winners_96_contig_list.txt", header = F)
contigs <- contigs$V1

contigs <- sub(pattern = ">", replacement = "", x = contigs)

genomes <- sub(pattern = "_scaffold.*$", replacement = "", x = contigs)

# ---- make coverM genome-reference ---- 

# genome.ref <- data.frame("genome" = genomes, "contig" = contigs)
# 
# head(genome.ref)
# tail(genome.ref)
# 
# write.table(x = genome.ref, file = "data/2023-01-25_coverM_files/input/coverM_genome_reference.txt", 
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- make inStrain stb ----

instrain.stb  <- data.frame("contig" = contigs, "genome" = genomes)

head(instrain.stb)
tail(instrain.stb)

write.table(x = instrain.stb, file = "data/2023-03-13_inStrain_on_drep96/input/instrain_stb.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ~
