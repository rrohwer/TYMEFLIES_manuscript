# RRR
# why they put it in a fucking ginormous excel file instead of .txt?
# why they make it basically impossible to find?
# That's just the JGI-way.


# Download the gold public dataset IDs from here: https://gold.jgi.doe.gov/downloads
# download limony IDs from SRA

library(readxl)
library(data.table)

file.1 <- "~/Dropbox/Research/2023-11-28_get_SRA_for_JGI_IDs/goldData.xlsx"
file.2 <- "~/Dropbox/Research/2023-11-28_get_SRA_for_JGI_IDs/sraBiomeImg.xlsx"

limony <- fread("~/Dropbox/Research/2023-11-28_get_SRA_for_JGI_IDs/limony_SraRunInfo.csv", colClasses = "character")

bigtyme <- fread("data/2023-12-08_combined_genome_info/sample_key.tsv", colClasses = "character")

# ---- fix a mistake in bigtyme ----

bigtyme[bigtyme$Extraction.Code == "rr 62", ]


# ---- get tymeflies existing accessions ----

foolsgold <- read_excel(path = file.1, sheet = "Sequencing Project", col_types = rep("text",23))

foolsgold <- as.data.table(foolsgold)

jgi.PID <- "504350"

gold <- foolsgold[`ITS PROPOSAL ID` %in% jgi.PID, .(`PROJECT GOLD ID`,`PROJECT NAME`,`SEQUENCING STRATEGY`,`ITS PROPOSAL ID`,`ITS SEQUENCING PROJECT ID`,`ITS SAMPLE ID`,`NCBI BIOPROJECT ACCESSION`,`NCBI BIOSAMPLE ACCESSION`,`SRA EXPERIMENT IDS`,`PROJECT STATUS`)]

gold <- gold[ ,.(`ITS SEQUENCING PROJECT ID`,`NCBI BIOPROJECT ACCESSION`,`NCBI BIOSAMPLE ACCESSION`,`SRA EXPERIMENT IDS`)]

bigtyme <- bigtyme[ ,`:=`(`NCBI Bioproject Accession` = NULL,`NCBI Biosample Accession` = NULL)]

key <- merge(x = bigtyme, y = gold, by.x = "ITS_SPID", by.y = "ITS SEQUENCING PROJECT ID")

key <- key[ ,.(sample,Extraction.Code,ITS_SPID,`NCBI BIOPROJECT ACCESSION`,`NCBI BIOSAMPLE ACCESSION`)]

colnames(key) <- c("TYMEFLIES.ID","limony.ID","JGI.ITS.SPID","Current.BioProject","Current.BioSample")

# lim.id <- key$limony.ID
# lim.num <- substr(x = lim.id, start = 4, stop = 100)
# index.na <- c(grep("tubes", x = lim.num, invert = F), grep("p", x = lim.num))
# 
# lim.num[nchar(lim.num) == 2] <- paste0("0", lim.num[nchar(lim.num) == 2])
# lim.num[nchar(lim.num) == 1] <- paste0("00", lim.num[nchar(lim.num) == 1])
# lim.num[index.na] <- NA
# 
# lim.let <- substr(x = lim.id, start = 1, stop = 2)
# lim.name <- paste0(lim.let, lim.num)
# lim.name[index.na] <- lim.id[index.na]
# 
# key$limony.ID <- lim.name

key[ ,limony.ID := sub(pattern = " ", replacement = "", x = limony.ID)]

# ---- get limony existing accessions ----

limony <- limony[ ,.(BioSample,SampleName)]
colnames(limony) <- c("Correct.BioSample","limony.ID")
limony[ ,limony.ID := sub(pattern = "p1$", replacement = "", x = limony.ID)]
limony[ ,limony.ID := sub(pattern = "p2$", replacement = "", x = limony.ID)]
limony[ ,limony.ID := sub(pattern = "p3$", replacement = "", x = limony.ID)]
limony[ ,limony.ID := sub(pattern = "p4$", replacement = "", x = limony.ID)]

# and what is the gen donor bioSample?

ugh <- readRDS(file = "pop/data/metadata/2020-07-10_data_entry_limony_and_TYMEFLIES.rds")
ugh <- ugh[ugh$sample.type == "generous.donor",]
# one is rr 998
ugh <- as.data.table(ugh)
ugh[ugh$Extraction.Code == "rr 998" & !is.na(ugh$Extraction.Code),]
# so rr998_p4
limony[limony.ID == "rr998"]
# should be SAMN29154387

# ---- match accessions ----

key <- merge(x = key, y = limony, by = "limony.ID", all.x = T, all.y = F)

key[grep(pattern = "gd", x = limony.ID), Correct.BioSample := "SAMN29154387"]

fwrite(x = key, file = "data/2024-02-09_fixing_SRA_numbers/Current_vs_Correct_BioSamples.tsv", sep = "\t")

# shit some of the limony ones need to be fixed too. for now just delete one and leave the other.