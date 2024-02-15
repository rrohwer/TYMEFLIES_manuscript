# RRR

# can use this website to get a map of metabolism pathways that are present:
# https://pathways.embl.de/ipath3.cgi?map=metabolic

# uploading the files manually is pretty time consuming, and you can't even export the animations!

# found this on programatically generating the figures: https://pathways.embl.de/tools.cgi
# and this on using POST with R: https://www.etiennebacher.com/posts/2023-05-09-making-post-requests-with-r/

# idea is to make a gif, where it lights up which genes are under selection in the different years

library(data.table)
library(lubridate)
library(httr) # get data from the web page
library(polite) # add delays when sending to web page to not do denial of service attack
library(magick) # convert svg to png, and make gifs from image files
library(ggplot2) # annotate the date onto the png. getting so hacky
library(ggpubr) # for background image

genome <- fread(file = "data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.ko-pathways.tsv")
outliers <- fread("data/2023-09-13_KEGG_annotations/example_files/ME2011-09-21_3300043464_group3_bin69.outlier_genes.tsv")
outliers <- unique(outliers$ko)
selection <- fread(file = "data/2023-08-05_example_gene_info_combined_MK_files/ME2011-09-21_3300043464_group3_bin69_gene_info_combined_MK.tsv.gz")
selection <- selection[MK.p.val <= .05 & mcdonald.kreitman < 1 & is.finite(mcdonald.kreitman), .(gene, sample, date, season)]
selection <- merge(x = selection, y = genome, by = "gene", all.x = TRUE, all.y = FALSE)

colnames(genome)
genome <- genome[ ,.(ko)]
genome[ ,width := "W10"]
genome[ ,color := "#F19DEE"]

fwrite(x = genome, file = "data/2023-09-27_metabolism_gif/base_genome.tsv", quote = F, col.names = F, sep = "\t")

genome.outliers <- copy(genome)
genome.outliers[ko %in% outliers, `:=`(color = "#0F16F5", width = "W30")]

# fwrite(x = genome.outliers, file = "data/2023-09-27_metabolism_gif/base_genome_with_outliers.tsv", quote = F, col.names = F, sep = "\t")

s = unique(selection$sample)[1]

for (s in unique(selection$sample)){
  pos <- selection[sample == s, ko]
  genome.day <- copy(genome)
  genome.day[ko %in% pos, `:=`(color = "#0F16F5", width = "W30")]
  # file.name <- sub(pattern = "IS_gene_info.*", replacement = "tsv", x = s)
  # fwrite(x = genome.day, file = file.path("data/2023-09-27_metabolism_gif",file.name), quote = F, col.names = F, sep = "\t", compress = "none")
  
  genome.day[ , color := sub(pattern = "#", replacement = "%23", x = color)]
  ko.table <- apply(X = genome.day, MARGIN = 1, FUN = paste, collapse = "%09")
  ko.table <- paste(ko.table, collapse = "%0D%0A")
  web.query <- paste0("selection=",ko.table,
                      "&default_opacity=1",
                      "&default_width=3",
                      "&default_radius=7",
                      "&default_color=%23aaaaaa",
                      "&background_color=%23ffffff",
                      "&tax_filter=",
                      "&map=metabolic",
                      "&export_type=svg",
                      "&export_dpi=120")
  dont.break.embl <- politely(POST, delay = 5) # post is the httr funciton: POST(url = "https://pathways.embl.de/mapping.cgi", body = web.query)
  web.response <- dont.break.embl(url = "https://pathways.embl.de/mapping.cgi", body = web.query)
  
  # connection <- file("~/Desktop/test1.svg", "wb") # This drops text for some reason?!?
  # writeBin(object = content(web.response), con = connection) 
  file.name <- sub(pattern = "IS_gene_info.*", replacement = "svg", x = s)
  write.table(x = rawToChar(content(web.response)), file = file.path("figures/2023-09-27_metabolism_gif/svgs",file.name), sep = "\n", quote = F, row.names = F, col.names = F)
  day.label <- sub("ME","",file.name)
  day.label <- sub("_.*$","",day.label)
  day.svg <- image_read_svg(path = file.path("figures/2023-09-27_metabolism_gif/svgs",file.name))
  day.png <- image_convert(image = day.svg, format = "png")
  day.png <- image_border(image = day.png, "white", "30x30")
  # image_annotate(image = day.png, text = day.label) # This doesn't work for some reason
  file.name <- sub(pattern = "IS_gene_info.*", replacement = "png", x = s)
  # image_write(image = day.png, path = file.path("figures/2023-09-27_metabolism_gif/pngs",file.name), format = "png")
  
  day.plot <- ggplot(data = data.frame("x" = 1:5,"y" = 1:5), aes(x = x, y = y)) + 
    theme_bw()+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = .5,size = 100)) + 
    background_image(raster.img = day.png) +
    labs(title = day.label)
  
  png(filename = file.path("figures/2023-09-27_metabolism_gif/pngs",file.name), width = 3774, height = 2250, units = "px") # orig svg aspect ratio: height='2250' width='3774'
  print(day.plot)
  dev.off()
}

# example query format:
# selection=K00010%09W10%09%23F19DEE%0D%0AK00027%09W10%09%23F19DEE%0D%0AK00057%09W10%09%23F19DEE&default_opacity=1&default_width=3&default_radius=7&default_color=%23aaaaaa&background_color=%23ffffff&tax_filter=&map=metabolic&
# web.query <- "selection=K00010+++W10+%23F19DEE&default_opacity=1&default_width=3&default_radius=7&default_color=%23aaaaaa&background_color=%23ffffff&tax_filter=&map=metabolic&export_type=svg&export_dpi=120"

# MAKe THE GIF TGIF

list.files(path="figures/2023-09-27_metabolism_gif/pngs/", pattern = '*.png', full.names = TRUE) %>% 
  image_read() %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=10) %>% # fps = frames per second
  image_write("figures/2023-09-27_metabolism_gif/acI-B_KOs_under_selection.gif", format = "gif") 

  
  