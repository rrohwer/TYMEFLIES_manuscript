# RRR


library(data.table)

gtdb <- fread(file = "data/2023-04-16_GTDB_accessions/1b-GTDB_download_Actinomycetia_Key.tsv")

ftp.paths <- gtdb[ ,ftp_path]

ftp.paths <- sub(pattern = "https", replacement = "ftp", x = ftp.paths)
ftp.paths <- sub(pattern = "$", replacement = "/\\*genomic.fna.gz", x = ftp.paths)

write.table(x = ftp.paths, file = "data/2023-04-16_GTDB_accessions/2_ftp_paths.txt", quote = F, row.names = F, col.names = F)
