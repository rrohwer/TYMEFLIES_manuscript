# RRR
# coverM needs to genome-definition file to pair with the reference file
# --reference 	PATH to FASTA file of contigs e.g. concatenated genomes
# --genome-definition		FILE containing list of genome_name<tab>contig lines to define the genome of each contig.

# I pulled out the list of contigs this way
# cat drep_winners_concat.fna | grep ">" > drep_winners_contig_list.txt
# read that in and parse it into the input file I need

contigs <- read.delim(file = "data/2023-02-01_coverM_files/input/drep_winners_contig_list.txt", header = F)
contigs <- contigs$V1

contigs <- sub(pattern = ">", replacement = "", x = contigs)

genomes <- sub(pattern = "_scaffold.*$", replacement = "", x = contigs)

genome.ref <- data.frame("genome" = genomes, "contig" = contigs)

instrain.stb  <- data.frame("contig" = contigs, "genome" = genomes)

head(genome.ref)
tail(genome.ref)

write.table(x = genome.ref, file = "data/2023-01-25_coverM_files/input/coverM_genome_reference.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(x = instrain.stb, file = "data/2023-02-01_coverM_files/input/instrain_stb.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
