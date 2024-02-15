# RRR- updated 4-4-23 with full inStrain results

library(lubridate)
library(ggplot2)
library(patchwork)
library(data.table)

# genome.info <- readRDS("data/2023-12-12_instrain_results_TIMEOUT/genome_info_combined.rds")
# bin.info <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")

genome.info <- fread(file = "data/2023-03-13_inStrain_on_drep96/processed_output/genome_info_combined.tsv.gz")
bin.info <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

bin.info <- as.data.table(bin.info)
bin.info <- bin.info[winner == T, ]

genome.info[ ,date := parse_date_time(x = as.character(date), orders = "ymd")]
genome.info[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]

colnames(genome.info)

# ---- CHOOSE 1 TO LOOK AT ----

plot.folder <- "figures/2023-02-16_instrain_first_looks/nucleotide_diversity/"
chosen.column <- "nucl_diversity"

plot.folder <- "figures/2023-02-16_instrain_first_looks/pop_ANI/"
chosen.column <- "popANI_reference"

plot.folder <- "figures/2023-02-16_instrain_first_looks/con_ANI/"
chosen.column <- "conANI_reference"

plot.folder <- "figures/2023-02-16_instrain_first_looks/linked_SNVs/"
chosen.column <- "linked_SNV_count"

# ---- 2 favorites ----

# 123	d__Bacteria	p__Actinobacteriota	c__Actinomycetia	o__Nanopelagicales	f__Nanopelagicaceae	g__AAA044-D11	s__
chosen.bin <- "ME2016-07-22_3300044628_group7_bin138"

# 106	d__Bacteria	p__Actinobacteriota	c__Actinomycetia	o__Nanopelagicales	f__Nanopelagicaceae	g__Nanopelagicus	s__Nanopelagicus sp000294575
chosen.bin <- "ME2011-09-21_3300043464_group3_bin69"

# 241	d__Bacteria	p__Cyanobacteria	c__Cyanobacteriia	o__Cyanobacteriales	f__Nostocaceae	g__Dolichospermum	s__Dolichospermum flosaquae
chosen.bin <- "ME2011-07-25_3300042383_group3_bin134"

my.info <- genome.info[genome == chosen.bin]
my.info <- my.info[order(date)]

# ---- look at all abundant ----

# subset top 200 for plots
med.cov.table <- dcast(data = genome.info, formula = genome ~ sample, value.var = "coverage_median")
med.cov.table <- as.matrix(x = med.cov.table, rownames = T)
overall.coverage <- apply(X = med.cov.table, MARGIN = 1, FUN = sum, na.rm = T)
overall.coverage <- overall.coverage[order(overall.coverage, decreasing = T)]
top.orgs <- overall.coverage[1:200]
barplot(overall.coverage)
abline(v = 200)

# ----

for (g in names(top.orgs)){
  my.info <- genome.info[genome == g]
  my.info$overall <- "Overall" # dummuy variable to make overall facet title
  
  taxonomy <- paste(bin.info[bin.full.name == g, c(phylum,class,order,family,genus,species,num.in.cluster)], collapse = " ") 
  short_plot_name <- gsub(pattern = "[pcofgs]__", replacement = "", x = paste(taxonomy,g))   # fucking powerpoint
  short_plot_name <- gsub(pattern = "33.*_group[12345678]_", replacement = "", x = short_plot_name)
  short_plot_name <- gsub(pattern = "ME", replacement = "", x = short_plot_name)
  
  # ----
  
  p.line.nuc.div <- ggplot(data = my.info, aes(x = date, y = .data[[chosen.column]]))+  
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(alpha = .5)+
    geom_point(aes(color = season), size = 1.5, alpha = .8)+
    scale_x_continuous(name = element_blank(), breaks = parse_date_time(seq.int(2000,2020,5),"y"), labels = seq.int(2000,2020,5))+
    guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5, alpha = 1)))
  
  p.line.cov.med <- ggplot(data = my.info, aes(x = date, y = coverage_median))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(alpha = .5)+
    geom_point(aes(color = season), size = 1.5, alpha = .8)+
    scale_x_continuous(name = element_blank(), breaks = parse_date_time(seq.int(2000,2020,5),"y"), labels = seq.int(2000,2020,5))+
    guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5, alpha = 1)))
  
  p.box.nuc.div <- ggplot(data = my.info, aes(x = season, y = .data[[chosen.column]]))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_boxplot(aes(fill = season), show.legend = F)+
    theme(axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank())
  
  p.box.cov.med <- ggplot(data = my.info, aes(x = season, y = coverage_median))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_boxplot(aes(fill = season), show.legend = F)+
    theme(axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank())

  p.box.inv.nucl.div <- ggplot(data = my.info, aes(x = invasion, y = .data[[chosen.column]]))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("steelblue","orange2","red3"))+
    scale_fill_manual(values = c("steelblue","orange2","red3"))+
    geom_boxplot(aes(fill = invasion), show.legend = F)+
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title.x = element_blank())+
    facet_wrap(~season, nrow = 1)
  
  p.box.inv.nucl.div.overall <- ggplot(data = my.info, aes(x = invasion, y = .data[[chosen.column]]))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("steelblue","orange2","red3"))+
    scale_fill_manual(values = c("steelblue","orange2","red3"))+
    geom_boxplot(aes(fill = invasion), show.legend = F)+
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title = element_blank())+
    facet_wrap(~overall)

  p.box.inv.cov.med <- ggplot(data = my.info, aes(x = invasion, y = coverage_median))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("steelblue","orange2","red3"))+
    scale_fill_manual(values = c("steelblue","orange2","red3"))+
    geom_boxplot(aes(fill = invasion), show.legend = F)+
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title.x = element_blank())+
    facet_wrap(~season, nrow = 1)
  
  p.box.inv.cov.med.overall <- ggplot(data = my.info, aes(x = invasion, y = coverage_median))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("steelblue","orange2","red3"))+
    scale_fill_manual(values = c("steelblue","orange2","red3"))+
    geom_boxplot(aes(fill = invasion), show.legend = F)+
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title = element_blank())+
    facet_wrap(~overall)

  p.line.nuc.div.year <- ggplot(data = my.info, aes(x = yday, y = .data[[chosen.column]]))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_smooth(color = "black", se = T, span = .4)+
    geom_point(aes(color = season, group = year), size = 3, alpha = .7, stroke = 0)+
    geom_smooth(color = "black", se = F, span = .4)+
    scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))+
    guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5, alpha = 1)))
  
  p.line.cov.med.year <- ggplot(data = my.info, aes(x = yday, y = coverage_median))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
    geom_line(aes(group = year), alpha = .2)+
    geom_smooth(color = "black", se = T, span = .4)+
    geom_point(aes(color = season, group = year), size = 3, alpha = .7, stroke = 0)+
    geom_smooth(color = "black", se = F, span = .4)+
    scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))+
    guides(color = guide_legend(title = "Season", override.aes = c(shape = 15, size = 5, alpha = 1)))
  
  p.nuc.inv.panel <- ggplot(data = my.info, aes(x = invasion, y = .data[[chosen.column]]))+ 
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
  
  p.cov.inv.panel <- ggplot(data = my.info, aes(x = invasion, y = coverage_median))+ 
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
  
  p.time.series <- p.line.cov.med + p.cov.inv.panel + p.line.nuc.div + p.nuc.inv.panel + plot_layout(guides = "collect", widths = c(4,1)) +
    plot_annotation(title = g, subtitle = taxonomy)
  
  p.year.overlay <- p.line.cov.med.year + p.box.cov.med + p.line.nuc.div.year + p.box.nuc.div + plot_layout(guides = "collect", widths = c(4,1)) +
    plot_annotation(title = g, subtitle = taxonomy)
  
  p.invasion.boxplots <- p.box.inv.cov.med + p.box.inv.cov.med.overall + p.box.inv.nucl.div + p.box.inv.nucl.div.overall + 
    plot_layout(nrow = 2, widths = c(6,1))+
    plot_annotation(title = g, subtitle = taxonomy)
  
  
  pdf(file = file.path(plot.folder,"timeline", paste(short_plot_name," - time_series.pdf")), width = 9.75, height = 6)
  print(p.time.series)
  dev.off()
  
  pdf(file = file.path(plot.folder, "time_overlay", paste(short_plot_name,"- year_overlay.pdf")), width = 9.75, height = 6)
  print(p.year.overlay)
  dev.off()
  
  pdf(file = file.path(plot.folder, "invasion_boxplots", paste(short_plot_name,"- inv_boxplots.pdf")), width = 9.75, height = 6)
  print(p.invasion.boxplots)
  dev.off()
  
  # ----
  
}







