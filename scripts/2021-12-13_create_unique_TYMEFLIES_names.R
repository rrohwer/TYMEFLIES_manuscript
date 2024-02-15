# RRR
# MY names are unique to the sample, and meaningful. 
# But will be more useful as filenames if they order by date
# the SPID is the sequencing project ID that is unique to each sample submission
# But when something is resequenced OR re-annotated it can have a duplicate SPID
# The taxon ID is unique to the assembly, so it is unique to the time it went through the pipeline
# So to keep my filenames unambiguous, 
# 1. search for them with SPID
# 2. manually choose which filtered reads to use when there are duplicates, so that it matches the taxon ID  
# 3. name them with the taxon ID that they correspond to later in the pipeline
# 4. ALSO name them with my own modified sample name:
#       ME16Jun2003D8 --> ME2003-06-16D8
# So filenames will look like this: ME2003-06-16D8_3300042865_filtered-reads.fastq.gz

# ---- set-up ----

library(readxl)
library(lubridate)

tf <- read_excel(path = "data/2021-12-13_IMG_IDs_for_504350/BIGTYME.xlsx", na = "NA")

lim <- readRDS(file = "~/Desktop/pop/data/limony-IRD/2021-08-25_processing/2C_key_samples_6000.rds")

ext <- read_excel("~/Dropbox/Research/Mendota_Extraction_Metadata/limony_and_TYMEFLIES/2018-04-01_choosing_filters_sorting_filters/2021-06-17_data_entry_limony_and_tymeflies.xlsx", na = "NA")

save.as <- "data/2021-12-13_IMG_IDs_for_504350/bigtyme_edits.csv"

# ---- fix extraction codes ----

tf$order <- 1:nrow(tf)

tf$Extraction.Code <- sub(pattern = "GENDONOR-", replacement = "gd", x = tf$Extraction.Code)

code.start <- substr(x = tf$Extraction.Code, start = 1, stop = 2)
code.end <- substr(x = tf$Extraction.Code, start = 3, stop = 10)
code.end <- as.numeric(code.end)
code.end <- as.character(code.end)
code.end[is.na(code.end)] <- sub(pattern = "l", replacement = "", x = substr(x = tf$Extraction.Code[is.na(code.end)], start = 3, stop = 10))

cbind(tf$Extraction.Code, paste(code.start, code.end))

tf$Extraction.Code <- paste(code.start, code.end)


# ---- get dates ----

tf$Sample.Name
tf$date <- substr(x = tf$Sample.Name, start = 3, stop = 11)
tf$date <- parse_date_time(tf$date, orders = "dmy")
tf$Year <- year(tf$date)
tf$Month <- month(tf$date)
tf$Day <- day(tf$date)

# ---- add filter codes and plate number ----

submitted <- data.frame("name" = ext$TYMEFLIES.name, "filt" = ext$Filter.Code, "ord" = ext$TYMEFLIES.order, "plate" = ext$TYMEFLIES.plate, "mystery" = ext$JGI.Sample.ID, "rep" = ext$Biological.Replicate)
submitted$name
submitted <- submitted[!is.na(submitted$name), ]
submitted <- submitted[!duplicated(submitted$name), ] # lines copied when 2 limony submissions of the sample



str(submitted)
submitted$name <- paste0("TYMEFLIES-",submitted$name)
cbind(submitted$name[1:5], tf$TYMEFLIES.JGI.Name[1:5])

submitted <- merge(x = submitted, y = tf, by.x = "name", by.y = "TYMEFLIES.JGI.Name", all = T)

submitted$name
submitted$ord

index.extras <- which(duplicated(submitted$ord))
index.extras <- submitted$ord %in% submitted$ord[index.extras]
# View(submitted[index.extras, ]) # ww0233 is double in IMG too
# now the dup has been deleted, so I can delete too in the excel file.

# which(is.na(submitted$ITS_SPID))
# submitted[180, ] # TYMEFLIES-ME12Jun2018-rr0366 not in IMG, added to spreadsheet...

tf <- tf[order(tf$order), ]
submitted <- submitted[order(submitted$order), ]

all.equal(submitted$name, tf$TYMEFLIES.JGI.Name) 

head(submitted)
tf$Filter.Code <- submitted$filt
tf$Submission.Plate <- submitted$plate
cbind(tf$`IMG Submission ID`, submitted$mystery)
tf$Submission.SPITS.ID <- submitted$mystery
tf$Submission.Order <- submitted$ord
head(tf)

tf$bio.rep <- submitted$rep

# ---- get tymeflies names ----

index <- grep(pattern = "S$", x = tf$Sample.Name, value = F)
tf$Sample.Name[index] <- sub(pattern = "S$", replacement = "D0", x = tf$Sample.Name[index])

tf$TYMEFLIES.name <- tf$Sample.Name

last.pf.date = "08Sep2003"
last.pf.date <- parse_date_time(x = last.pf.date, orders = "dmy")
index <- order(tf$date)
tf <- tf[index, ]

index.old <- which(tf$date <= last.pf.date)
index.ww <- which(tf$bio.rep == "W") # they all are distinguished with a W that I added manually in data entry, all uppercase
index.pf <- setdiff(x = index.old, y = index.ww)

tf$TYMEFLIES.name[index.pf] <- paste0(tf$TYMEFLIES.name[index.pf], "pf")

tf$TYMEFLIES.name

lim$in.R.colnames

# format so filenames go in order by date:

my.date <- substr(x = tf$TYMEFLIES.name, start = 3, stop = 11)
my.date <- parse_date_time(x = my.date, orders = "dmy", tz = "Etc/GMT-5")
my.notes <- substr(x = tf$TYMEFLIES.name, start = 12, stop = 100)
my.format <- stamp(x = c("2000-01-01","2002-12-31","1999-06-15"), orders = "ymd")
my.date <- my.format(my.date)
new.names <- paste0("ME",my.date, my.notes)
tf$TYMEFLIES.name <- new.names

# ---- export csv, paste in edited columns to formatted excel big file ----

tf <- tf[order(tf$order), ]

write.csv(x = tf, file = save.as, row.names = F, quote = F)
cat(save.as)
