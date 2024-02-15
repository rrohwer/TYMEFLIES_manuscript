# RRR
# some pretty examples of microdiversity change over time
# use abundant acI-B as examples: 

# Notes:
# Ca. Planktophila dulcis 	acI-A1		
# Ca. Planktophila sulfonica 	acI-A1
# Ca. Planktophila versatilis 	acI-A1
# Ca. Planktophila lacus 	acI-A1
# Ca. Planktophila rubra 	acI-A1
# Ca. Planktophila limnetica	acI-A2		
# Ca. Planktophila aquatilis 	acI-A4
# Ca. Planktophila vernalis 	acI-A7

# Ca. Planktophila limnetica 	Phila

# Ca. Nanopelagicus limnes 	acI-B1
# Ca. Nanopelagicus hibericus 	acI-B1
# Ca. Nanopelagicus abundans 	acI-B1

library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

genome.info <- readRDS("data/2023-12-12_instrain_results_TIMEOUT/genome_info_combined.rds")
bin.info <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")
bin.info <- bin.info[bin.info$winner, ]

# ----

g <- "ME2015-06-27_3300042558_group6_bin43"
my.info <- genome.info[genome.info$genome == g, ]
my.info$year <- year(my.info$date)
my.info$yday <- yday(my.info$date)
my.info$overall <- "Overall" # dummuy variable to make overall facet title

taxonomy <- paste(bin.info[bin.info$bin.full.name == g, c("phylum","class","order","family","genus","species","num.in.cluster")], collapse = " ") 
short_plot_name <- gsub(pattern = "[pcofgs]__", replacement = "", x = paste(taxonomy,g))   # fucking powerpoint
short_plot_name <- gsub(pattern = "33.*_group[12345678]_", replacement = "", x = short_plot_name)
short_plot_name <- gsub(pattern = "ME", replacement = "", x = short_plot_name)

p.nuc <- ggplot(data = my.info, aes(x = yday, y = nucl_diversity))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = 0)+
  geom_smooth(color = "black", se = T)+
  geom_point(aes(color = season, group = year), size = 3, alpha = .8)+
  geom_smooth(color = "black", se = F)+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))+
  scale_y_continuous(name = "Microdiversity\n(nucleotide diversity)")+
  guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5)))

p.cov <- ggplot(data = my.info, aes(x = yday, y = coverage_median))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = 0)+
  geom_smooth(color = "black", se = T)+
  geom_point(aes(color = season, group = year), size = 3, alpha = .8)+
  geom_smooth(color = "black", se = F)+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))+
  scale_y_continuous(name = "Abundance\n(median coverage)")+
  guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5)))

p.nuc.box <- ggplot(data = my.info, aes(x = season, y = nucl_diversity))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
  geom_boxplot(aes(fill = season), show.legend = F)+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank())

p.cov.box <- ggplot(data = my.info, aes(x = season, y = coverage_median))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
  geom_boxplot(aes(fill = season), show.legend = F)+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank())


pdf(file = paste0("figures/2023-02-19_for_IGB_talk/acI-B_year-overlay_",g,".pdf"), width = 9, height = 4)

p.cov + p.cov.box + p.nuc + p.nuc.box + plot_layout(guides = "collect", widths = c(4,1)) +
  plot_annotation(title = "Actinobacteria Nanopelagicus (acI-B bin 43)", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()

# ----

g <- "ME2011-09-21_3300043464_group3_bin69"

my.info <- genome.info[genome.info$genome == g, ]
my.info$year <- year(my.info$date)
my.info$yday <- yday(my.info$date)
my.info$overall <- "Overall" # dummuy variable to make overall facet title

taxonomy <- paste(bin.info[bin.info$bin.full.name == g, c("phylum","class","order","family","genus","species","num.in.cluster")], collapse = " ") 
short_plot_name <- gsub(pattern = "[pcofgs]__", replacement = "", x = paste(taxonomy,g))   # fucking powerpoint
short_plot_name <- gsub(pattern = "33.*_group[12345678]_", replacement = "", x = short_plot_name)
short_plot_name <- gsub(pattern = "ME", replacement = "", x = short_plot_name)

p.nuc <- ggplot(data = my.info, aes(x = date, y = nucl_diversity))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(alpha = 1)+
  geom_point(aes(color = season), size = 1.5, alpha = .8)+
  scale_x_continuous(name = element_blank(), breaks = parse_date_time(seq.int(2000,2020,5),"y"), labels = seq.int(2000,2020,5))+
  scale_y_continuous(name = "Microdiversity\n(nucleotide diversity)")+
  guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5)))

p.cov <- ggplot(data = my.info, aes(x = date, y = coverage_median))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(alpha = 1)+
  geom_point(aes(color = season), size = 1.5, alpha = .8)+
  scale_x_continuous(name = element_blank(), breaks = parse_date_time(seq.int(2000,2020,5),"y"), labels = seq.int(2000,2020,5))+
  scale_y_continuous(name = "Abundance\n(median coverage)")+
  guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5)))

p.nuc.inv <- ggplot(data = my.info, aes(x = invasion, y = nucl_diversity))+ 
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("steelblue","orange2","red3"),labels = c("None","Spiny Water Flea","Zebra Mussel"))+
  scale_color_manual(values = c("steelblue","orange2","red3"),labels = c("None","Spiny Water Flea","Zebra Mussel"))+
  geom_boxplot(aes(fill = invasion), show.legend = F)+
  geom_point(aes(color = invasion), show.legend = T, alpha = 0)+ # for legend
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank())+
  guides(color = guide_legend(title = "Invasion", override.aes = c(shape = 15, alpha = 1, size = 5)))+
  guides(fill = FALSE)

p.cov.inv <- ggplot(data = my.info, aes(x = invasion, y = coverage_median))+ 
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("steelblue","orange2","red3"),labels = c("None","Spiny Water Flea","Zebra Mussel"))+
  scale_color_manual(values = c("steelblue","orange2","red3"),labels = c("None","Spiny Water Flea","Zebra Mussel"))+
  geom_boxplot(aes(fill = invasion), show.legend = F)+
  geom_point(aes(color = invasion), show.legend = T, alpha = 0)+ # for legend
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank())+
  guides(color = guide_legend(title = "Invasion", override.aes = c(shape = 15, alpha = 1, size = 5)))+
  guides(fill = FALSE)


pdf(file = paste0("figures/2023-02-19_for_IGB_talk/acI-B_invasion-change_",g,".pdf"), width = 9, height = 4)

p.cov + p.cov.inv + p.nuc + p.nuc.inv + plot_layout(guides = "collect", widths = c(4,1)) +
  plot_annotation(title = "Actinobacteria Nanopelagicus (acI-B bin 69)", 
                  theme = theme(plot.title = element_text(hjust = .5)))
dev.off()

