# RRR
# change from bin.1.fa to ME2000-03-15pf_3300044539_group1_bin1.fna

# ---- import -----

groups <- readRDS(file = "data/2022-06-27_binning_groups/mapping_groups_list.rds")
groups <- groups[-1]

# output of ls -R metabat2_bins with unzipped bin folders:
# $ pwd
# /work/08922/rrohwer/ls6/metabat2_bins
# $ ls -R > ../binlist.txt

bins <- read.delim(file = "data/2022-09-14_list_of_all_bins/binlist.txt")
bins <- bins[ ,1,drop = T]


# ---- functions ----

make.empty.list.structure <- function(ListNames){
  # the ListNames can be something like c("OTU", "kingdom","phylum","class","order","family/lineage","genus/clade","species/tribe")
  empty.list <- list(NULL)
  for (e in 1:length(ListNames)){
    empty.list[[e]] <- 0
    names(empty.list)[e] <- ListNames[e]
  }
  return(empty.list)
}

# ---- format bins ----

index.folders <- grep(x = bins, pattern = "^ME", value = F)
bins <- bins[-index.folders]

index.folders <- grep(x = bins, pattern = "^\\.\\/ME", value = F)
names.folders <- grep(x = bins, pattern = "^\\.\\/ME", value = T)
names.folders <- sub(pattern = "^\\.\\/", replacement = "", x = names.folders)
names.folders <- sub(pattern = ":$", replacement = "", x = names.folders)

bin.list <- make.empty.list.structure(ListNames = names.folders)

length(bin.list)
length(names.folders)
length(index.folders)

for (i in 1:length(bin.list)){
  start <- index.folders[i] + 1
  if (i == length(bin.list)){
    end <- length(bins)
  }else{
    end <- index.folders[i + 1] - 1
  }
  bin.list[[i]] <- bins[start:end]
}

# ---- parse into tables ----

all.bins <- unlist(x = bin.list, use.names = F) # using names it adds a count at the end of each name

all.bins <- NULL
all.names <- NULL
for (e in 1:length(bin.list)){
  num.bins <- length(bin.list[[e]])
  all.bins <- c(all.bins, bin.list[[e]])
  all.names <- c(all.names, rep(x = names(bin.list)[e], times = num.bins))
}

bin.table <- data.frame("sample" = all.names, "bin" = all.bins)
head(bin.table)

all.samples <- NULL
all.groups <- NULL
for (e in 1:length(groups)){
  num.samples <- length(groups[[e]])
  all.samples <- c(all.samples, groups[[e]])
  all.groups <- c(all.groups, rep(x = names(groups)[e], times = num.samples))
}

group.table <- data.frame("group" = all.groups, "sample" = all.samples)
head(group.table)

bin.table <- merge(x = bin.table, y = group.table, by = "sample", all = TRUE)
head(bin.table)

rename.table <- bin.table
rename.table$bin.num <- sub(pattern = "bin\\.", replacement = "", x = rename.table$bin)
rename.table$bin.num <- sub(pattern = "\\.fa", replacement = "", x = rename.table$bin.num)
rename.table$new.name <- paste0(rename.table$sample, "_", rename.table$group, "_bin", rename.table$bin.num,".fna")

bin.table <- rename.table[ ,-4]
colnames(bin.table)[4] <- "name"
head(bin.table)

# ---- export ----

saveRDS(object = bin.list, file = "data/2022-09-14_list_of_all_bins/binlist.rds")

saveRDS(object = bin.table, file = "data/2022-09-14_list_of_all_bins/bin_key.rds")

# ---- end ----
