# RRR

# need to rename the scaffolds from each bin to add the bin names to them.
# the bins were already re-named
# use the un-renamed example bin file to test though:
# userprefs <- "data/2022-08-02_example_bins"
# userprefs[2] <- "data/2022-10-18_test_scaffold_renaming"

# idev -t 2:00:00
# module load Rstats
# Rscript rename_scaffolds.R metabat2_bins metabat2_bins_renamed
# had to run as submitted script- took 02:57:36
# running submit script with more flexible regex tool 04:31:04



userprefs <- commandArgs(trailingOnly = TRUE)

metabat.folder <- userprefs[1]
new.metabat.folder <- userprefs[2]

my.samples <- list.dirs(path = metabat.folder, recursive = F, full.names = F)

for (s in my.samples){
  cat("working on",s,"\n")
  bin.folder <- s
  new.bin.folder <- file.path(new.metabat.folder,bin.folder)
  dir.create(new.bin.folder)
  
  all.bins <- list.files(path = file.path(metabat.folder,bin.folder))
  
  for (b in 1:length(all.bins)){
    my.bin <- all.bins[b]
    my.seqs <- read.table(file = file.path(metabat.folder, bin.folder, my.bin))
    my.seqs <- my.seqs[ ,1,drop = T]
    my.seqs[1:5]
    index.scaffolds <- grep(pattern = ">",x = my.seqs, value = F)
    old.scaffolds <-  grep(pattern = ">",x = my.seqs, value = T)
    old.scaffolds <- sub(pattern = ">",replacement = "",x = old.scaffolds)
    name.addition <- sub(pattern = "\\.[[:alpha:]]*$", replacement = "", x = my.bin)
    new.scaffolds <- paste0(">", name.addition, "_", old.scaffolds)
    cbind(my.seqs[index.scaffolds],new.scaffolds)
    my.seqs[index.scaffolds] <- new.scaffolds
    write.table(x = my.seqs, file = file.path(new.bin.folder, my.bin), quote = F, row.names = F, col.names = F)
  }
}





