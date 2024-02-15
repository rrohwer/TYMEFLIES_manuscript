# RRR

# samtools doesn't use all the threads
# but limited by memory use- max 100 GB?
# 250 GB per node, so can do 2 at once
# but also, keep within the mapping groups

mapping.groups <- readRDS("data/2022-06-27_binning_groups/mapping_groups_list.rds")
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

# ---- get samtools groups ----
samtools.groups <- make.empty.list.structure(ListNames = names(mapping.groups))
for (m in names(mapping.groups)){
  my.map <- mapping.groups[[m]]
  my.sam <- make.empty.list.structure(ListNames = paste0("group",1:ceiling(length(my.map) / 2)))
  s <- 1
  for (i in seq.int(from = 1, to = length(my.map), by = 2)){
    my.sam[[s]] <- my.map[c(i, i +1)]
    if (any(is.na(my.sam[[s]]))){
      my.sam[[s]] <- my.sam[[s]][1]
    }
    s <- s + 1
  }
  samtools.groups[[m]] <- my.sam
}

# ---- get submission checklist ----

# wait I don't need to do this at all, because can run all per assembly at once and just specify different node breakdown
# it is all per_assembly not per group. the groups just determine the cross mapping, but each assembly is submitted separately.