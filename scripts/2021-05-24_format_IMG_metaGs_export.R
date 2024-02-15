# RRR
# I downloaded the metaG (genome table) stats from IMG, so take big picture looks at them
# I know their quality, size, # bins, etc
# Note on IDs:
#  "ITS SP ID" same as "JGI Project ID"
#  "taxon_oid" same as "IMG Genome ID"

# ---- set-up and import ----

library(readxl)
library(tidyr)

clean.rds.output <- "data/2021-05-24_metaG_stats/all_img_assemblies.rds"
clean.tsv.output <- "data/2021-05-24_metaG_stats/all_img_assemblies.tsv"

ids.only.output.its.sp.id <- "data/2021-05-24_metaG_stats/ITS_SP_IDs.txt"
ids.only.output.taxon.id <- "data/2021-05-24_metaG_stats/Taxon_IDs.txt"

metaGs <- read_excel(path = "data/2021-05-24_metaG_stats/2021-05-21_all_sequenced_504350.xlsx", 
                     col_types = c("text","text","text","text","text","text","text","text","text","text","text","text","text","numeric","numeric","numeric","numeric"))
colnames(metaGs)
my.spreadsheet <- readRDS(file = "data/2021-05-24_metaG_stats/2020-07-10_data_entry_limony_and_TYMEFLIES.rds")
colnames(my.spreadsheet)

# ---- subset my spreadsheet of limony/tymeflies to only submitted tymeflies samples ----

index1 <- which(!is.na(my.spreadsheet$JGI.Sample.ID))
index2 <- which(!is.na(my.spreadsheet$TYMEFLIES.name))
all.equal(index1,index2)

my.spreadsheet <- my.spreadsheet[index1, ]

colnames(my.spreadsheet)

#
# ---- look at img table- wtf is up with these duplicate columns with different names?? ----

all.equal(metaGs[ ,1], metaGs[ ,7])
which(names(metaGs[ ,1]) != names(metaGs[ ,7]))
# OK these are the same
names(metaGs[1,1]) # "taxon_oid"
names(metaGs[1,7]) # "IMG Genome ID"

all.equal(metaGs[ ,13], metaGs[ ,9])
which(names(metaGs[ ,13]) != names(metaGs[ ,9]))
names(metaGs[1,13]) # "ITS SP ID"
names(metaGs[1,9]) # "JGI Project ID / ITS SP ID"
index <- which(metaGs[ ,13] != metaGs[ ,9])
# View(metaGs[index, c(9,13)])
# OK, all the diffs are the 9 is zero. why??
# whatever go with 13

# remove confusing duplicate columns:
metaGs <- metaGs[ ,-c(9,7)]

# and take all these spaces and weird characters out of the names!
colnames(metaGs)
colnames(metaGs) <- gsub(pattern = " ", replacement = ".", x = colnames(metaGs)) %>%
  gsub(pattern = "\\*", replacement = "", x = .) %>%
  gsub(pattern = "/", replacement = "", x = .) %>%
  gsub(pattern = "\\.\\.\\.\\.", replacement = ".", x = .)
colnames(metaGs)


# # ---- merge using the submission IDs ... or don't ----
# 
# colnames(my.spreadsheet)
# colnames(metaGs) <- gsub(pattern = " ", replacement = ".", x = colnames(metaGs))
# colnames(metaGs)
# 
# my.spreadsheet$JGI.Sample.ID[1:5]
# metaGs$IMG.Submission.ID[1:5]
# length(unique(c(my.spreadsheet$JGI.Sample.ID,metaGs$IMG.Submission.ID))) #933 these don't match wtf
# str(my.spreadsheet$JGI.Sample.ID)
# str(metaGs$IMG.Submission.ID)
# cbind(sort(my.spreadsheet$JGI.Sample.ID), sort(metaGs$IMG.Submission.ID))
# 
# giant.table <- merge(x = my.spreadsheet, y = metaGs, by.x = "JGI.Sample.ID", by.y = "IMG.Submission.ID")
# 
# # OK I have no clue what fucking ID is in my spreadsheet, it's from the submission.
# # it doesn't match anything here and I just clicked alllll the fucking ids when I downloaded from img
# # god their IDs suck and a table of what they all are "doesn't

# ---- parse out my sample IDs from IMG table ----

metaGs$Genome.Name..Sample.Name

my.IDs <- metaGs[ ,5, drop =T] %>%
  sub(pattern = "Freshwater microbial communities from Lake Mendota, Madison, Wisconsin, United States - TYMEFLIES-", replacement = "", x = .) %>%
  sub(pattern = "Freshwater microbial communities from Lake Mendota, Wisconsin, United States - TYMEFLIES-", replacement = "", x = .)
my.IDs

metaGs <- cbind("TYMEFLIES.name" = my.IDs, metaGs)

# ---- fix names that don't match (the gen donors) ----

setdiff(x = my.spreadsheet$TYMEFLIES.name, y = metaGs$TYMEFLIES.name) # in my spreadsheet and not IMG
setdiff(x = metaGs$TYMEFLIES.name, y = my.spreadsheet$TYMEFLIES.name) # in IMG and not my spreadsheet
# rename the controls in in both- to match and to include a sample name now.
i <- which(my.spreadsheet$TYMEFLIES.name == "CONTROL-GENDONOR-p1")
my.spreadsheet$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0001"
i <- which(metaGs$TYMEFLIES.name == "CONTROL-GENDONOR")
metaGs$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0001"

i <- which(my.spreadsheet$TYMEFLIES.name == "CONTROL-GENDONOR-p2")
my.spreadsheet$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0002"
i <- which(metaGs$TYMEFLIES.name == "CONTROL-GENDONOR pl2")
metaGs$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0002"

i <- which(my.spreadsheet$TYMEFLIES.name == "CONTROL-GENDONOR-p3")
my.spreadsheet$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0003"
i <- which(metaGs$TYMEFLIES.name == "CONTROL-GENDONOR pl3")
metaGs$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0003"

i <- which(my.spreadsheet$TYMEFLIES.name == "CONTROL-GENDONOR-p4")
my.spreadsheet$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0004"
i <- which(metaGs$TYMEFLIES.name == "CONTROL-GENDONOR pl4")
metaGs$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0004"

i <- which(my.spreadsheet$TYMEFLIES.name == "CONTROL-GENDONOR-tb")
my.spreadsheet$TYMEFLIES.name[i] <- "ME08Nov2018GD-tb0001"
i <- which(metaGs$TYMEFLIES.name == "CONTROL-GENDONOR tubes")
metaGs$TYMEFLIES.name[i] <- "ME08Nov2018GD-tb0001"

i <- which(my.spreadsheet$TYMEFLIES.name == "CONTROL-GENDONOR-p5")
my.spreadsheet$TYMEFLIES.name[i] <- "ME08Nov2018GD-pl0005"

setdiff(x = my.spreadsheet$TYMEFLIES.name, y = metaGs$TYMEFLIES.name) # these ones were sent but didn't make it through the pipeline
setdiff(x = metaGs$TYMEFLIES.name, y = my.spreadsheet$TYMEFLIES.name) # good should be zero


# ---- merge and clean up ----

giant.table <- merge(x = my.spreadsheet, y = metaGs, by = "TYMEFLIES.name", all.x = T, all.y = T)

cbind(1:ncol(giant.table), colnames(giant.table))

assemblies <- data.frame("IMG.Taxon.ID" = giant.table$taxon_oid,
                         "ITS.SP.ID" = giant.table$ITS.SP.ID,
                         "IMG.Submission.ID" = giant.table$IMG.Submission.ID,
                         "ITS.Proposal.ID" = giant.table$ITS.Proposal.ID,
                         "GOLD.Analysis.Project.ID" = giant.table$GOLD.Analysis.Project.ID,
                         "GOLD.Sequencing.Project.ID" = giant.table$GOLD.Sequencing.Project.ID,
                         "Sequencing.Plate" = giant.table$TYMEFLIES.plate,
                         "TYMEFLIES.Name" = giant.table$TYMEFLIES.name,
                         "Sample.Name" = giant.table$Sample.Name,
                         "Filter.Code" = giant.table$Filter.Code,
                         "Extraction.Code" = giant.table$Extraction.Code,
                         "Assembly.Size" = giant.table$Genome.Size.assembled,
                         "Num.Scaffolds" = giant.table$Scaffold.Count.assembled,
                         "Num.Genes" = giant.table$Gene.Count.assembled,
                         "Num.Bins" = giant.table$Gene.Count.assembled,
                         "DNA.Extraction.ng.mL" = giant.table$DNA.ng.uL,
                         "Filter.Notes" = giant.table$Filter.Notes,
                         "Extraction.Notes" = giant.table$Extraction.Notes)
str(assemblies)

sp.ids.only <- assemblies$ITS.SP.ID # "ITS SP ID" same as "JGI Project ID"
taxon.ids.only <- assemblies$IMG.Taxon.ID # "taxon_oid" same as "IMG Genome ID"

# ---- export ----

# saveRDS(object = assemblies, file = clean.rds.output)
# write.table(x = assemblies, file = clean.tsv.output, quote = F, row.names = F)
# 
# write.table(x = sp.ids.only, file = ids.only.output.its.sp.id, quote = F, row.names = F)
# write.table(x = taxon.ids.only, file = ids.only.output.taxon.id, quote = F, row.names = F)


# ~end~

# unnecessary parsing:
# Sample.Name <- sub(pattern = "-.*$", replacement = "", x = my.IDs)
# Sample.Name
# TF.ex.Code <- sub(pattern = "^.*-", replacement = "", x = my.IDs)
# rr <- substr(x = TF.ex.Code, start = 1, stop = 2)
# ex.num <- substr(x = TF.ex.Code, start = 3, stop = 7) %>%
#   as.numeric(.) # these NAs are gen donors, take care of in next code block
# Extraction.Code <- paste(rr, ex.num)
# 
# grep(pattern = "ME", x = Sample.Name, ignore.case = T, value = T, invert = T) # p5 control is missing
# index <- which(Sample.Name == "CONTROL")
# Sample.Name[index] <- "ME08Nov2018GD"
# Extraction.Code[index] <- TF.ex.Code[index] %>%
#   sub(pattern = "^GENDONOR", replacement = "gd", x = .)

