# RRR

# There is no master table of K numbers and their associated pathways
# Each K number can be associated with multiple pathways (b/c one enzyme can be involved in multiple pathways)
# But there is a website where you can ENTER the K number and LOOKUP all of it's pathways. here:

# https://www.kegg.jp/kegg/rest/keggapi.html#list

# And on that same website, you can download the pathway code and its description as a table.
# (That is the pathways.txt file in the KEGG_lists folder)

# So here's the plan:
# Run the annotation of all genomes with Kofamscan
# Get a vector of unique K numbers
# loop through all those K numbers with the API interface to get their pathways (~slow~)
#   do this using the httr2 package. example syntax for looking up K03889 is:
#   x <- request(base_url = "https://rest.kegg.jp/link/pathway/K03889") # get the request
#   y <- req_perform(stupid)                                            # make the request
#   z <- resp_body_string(resp = stupider)                              # interp the response as a string
#   my.tab <- fread(text = z, sep = "\t", header = F)                   # parse the string into a table
# Run keggdecoder on each genome to get a list of pathways present in it
# For each gene with a K number, assign the pathway that ALSO exists in its genome



library(httr2)
library(data.table)

pathways <- fread(file = "data/2023-09-13_KEGG_annotations/KEGG_lists/pathway_list.txt", sep = "\t", header = F)

kofamscan.example <- fread(file = "data/2023-09-13_KEGG_annotations/ME2011-09-21_3300043464_group3_bin69.kofamscan.tsv", sep = "\t", header = T, )
kofamscan.example <- kofamscan.example[-1, ]
colnames(kofamscan.example) <- c("is.signif","gene","KO","threshold","score","e.value","KO.definition")

unique.KOs <- unique(kofamscan.example$KO)
# for testing
unique.KOs <- unique.KOs[1:10]
unique.KOs <- paste0("https://rest.kegg.jp/link/pathway/",unique.KOs)

key.list <- list(NULL)
for (k in unique.KOs){
  my.tab <- request(base_url = k) %>%
    req_perform() %>%
    resp_body_string() %>%
    fread(text = ., sep = "\t", header = F)
 my.tab <- my.tab[grep(x = V2, pattern = "path:ko", value = F, invert = T)] 
 my.tab[ ,`:=`(V1 = sub(pattern = "ko:", replacement = "",x = V1), V2 = sub(pattern = "path:", replacement = "", x = V2))]
 key.list <- c(key.list, my.tab)
}

# OH MY GOD!!!  I CAN JUST DOWNLOAD THE WHOLE FUCKING TABLE THIS WAY:

# https://rest.kegg.jp/link/pathway/ko

# no need for the https lookups. OK moving on.
