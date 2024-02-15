# RRR
# one bin at a time right now.

# ---- set up ----

library(seqinr)

input <- commandArgs(trailingOnly = TRUE)
binid <- input[1]
foldergff <- input[2]
folderfna <- input[3]
foldernew <- input[4]

# test
cat("\ncomment out test files\n")
binid <- "3300042290_70"
foldergff <- "finding-data/2022-01-05_actino_bins"
folderfna <- "finding-data/2022-01-05_actino_bins"
foldernew <- "finding-data/2022-01-05_actino_bins"

filegff <- file.path(foldergff,paste0(binid,".gff"))
filefna <- file.path(folderfna,paste0(binid,".fna"))
filenew <- file.path(foldernew,paste0(binid,"_16S.fasta"))

gff <- read.delim(file = filegff, header = F, sep = "\t")
fna <- read.fasta(file = filefna, seqtype = "DNA", as.string = T, forceDNAtolower = F)

# ---- do things ----

fna <- data.frame("SeqID" = names(fna), "seq" = unlist(fna))
colnames(gff) <- c("SeqID","Version","FeatureType","Start","End","Score", "Strand","Phase","Details")

index.rRNA <- grep(pattern = "rRNA",x = gff$FeatureType, ignore.case = T, value = F)
index.16S <- grep(pattern = "16S",x = gff$Details, ignore.case = T, value = F)
index.16SrRNA <- intersect(index.rRNA, index.16S)

gff <- gff[index.16SrRNA, ]

fna <- merge(x = fna, y = gff, by = "SeqID", all.x = F, all.y = T)

fna$seq <- substr(x = fna$seq, start = fna$Start, stop = fna$End)

fna$SeqID <- paste(binid, fna$SeqID, 1:length(fna$SeqID), sep = "-")
fna$SeqID <- paste0(">",fna$SeqID)
fna <- fna[ ,1:2]

# ----
cat("Making file: ", filenew, "\n")
write.table(x = fna, file = filenew, quote = F, sep = "\n", row.names = F, col.names = F)
