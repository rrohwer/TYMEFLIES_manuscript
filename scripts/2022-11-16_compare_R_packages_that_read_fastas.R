
library(seqinr)

x <- read.fasta(file = "data/2022-08-02_example_bins/ME2000-03-15pf_3300044539/bin.1.fa")
# does each scaffold as a list, each base as a vector element
# seems inefficient

x <- read.fasta(file = "data/2022-08-02_example_bins/ME2000-03-15pf_3300044539/bin.1.fa", as.string = T)
x
# well that's a little better, sequence is a single string

library(phylotools)

x <- read.fasta(file = "data/2022-08-02_example_bins/ME2000-03-15pf_3300044539/bin.1.fa", clean_name = F)
str(x)
# this is a better table-format

nchar(x$seq.text)


library(Biostrings)
# way harder to install and tons of dependencies (so phylotools seems lighter)

x <- readDNAStringSet(filepath = "data/2022-08-02_example_bins/ME2000-03-15pf_3300044539/bin.1.fa", format = "fasta")
x
# oh wow it is like a special colored format
str(x)
# with like all this other info


# SO, for simply reading in sequences, phylotools seems best
# but if you want to manipulate the sequence then Biostrings seems to have a lot of functionality
# I'm just not impressed by the list-structure choice of seqinr, but makes sense for having each character a separate element