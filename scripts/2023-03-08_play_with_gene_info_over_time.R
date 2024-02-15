# RRR

library(data.table)
library(lubridate)
library(ggplot2)
library(patchwork)

# fread and fwrite would be faster for large files, and can still zip IF Open MP enabled!
genes <- fread(file = "data/2023-02-12_instrain_results_TIMEOUT/gene_info_combined.tsv.gz")

tax <- readRDS(file = "data/2022-11-15_bin_stats/all_bins.rds")

genes[ ,date := parse_date_time(x = as.character(date), orders = "ymd", tz = "Etc/GMT-5")]
genes[ ,season := factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall"))]
colnames(genes)

my.genome <- unique(genes$genome)[1]
my.genome <- "ME2015-06-27_3300042558_group6_bin43"

one.genome <- genes[genome == my.genome]
colnames(one.genome)

ggplot(data = one.genome, aes(x = date, y = divergent_site_count))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_line(aes(group = gene), alpha = .2)+
  scale_x_continuous(breaks = parse_date_time(x = 2000:2020, orders = "y"), labels = 2000:2020)+
  geom_vline(xintercept = parse_date_time(x = 2000:2020, orders = "y"), color = "yellow")

p.cov <- ggplot(data = one.genome, aes(x = yday, y = coverage))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = adjustcolor(c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"),alpha.f = .2))+
  geom_point(aes(color = season), alpha = .2)+
  geom_smooth(aes(group = gene), se = F, col = adjustcolor("black",alpha.f = .2))+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))

p.nuc.div <- ggplot(data = one.genome[!is.na(nucl_diversity) & nucl_diversity != 0  & nucl_diversity < .2], aes(x = yday, y = nucl_diversity))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = adjustcolor(c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"),alpha.f = .2))+
  geom_point(aes(color = season), alpha = .2)+
  # geom_smooth(aes(group = gene), se = F, col = adjustcolor("black",alpha.f = .2), span = 1)+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))
  # coord_trans(y = "log")

p.cov + p.nuc.div + plot_layout(nrow = 2, guides = "collect")


one.genome.average <- one.genome[ ,.(nuc.div = mean(nucl_diversity, na.rm = T), cov = mean(coverage, na.rm = T)), by = c("year","yday","season")]

p.cov <- ggplot(data = one.genome.average, aes(x = yday, y = cov))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = adjustcolor(c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"),alpha.f = 1))+
  geom_point(aes(color = season), alpha = 1, size = 3)+
  geom_smooth(se = F, col = adjustcolor("black",alpha.f = .5))+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))

p.nuc.div <- ggplot(data = one.genome.average[!is.na(nuc.div) & nuc.div != 0  & nuc.div < .2], aes(x = yday, y = nuc.div))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = adjustcolor(c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"),alpha.f = 1))+
  geom_point(aes(color = season), alpha = 1, size = 3)+
  geom_smooth(se = F, col = adjustcolor("black",alpha.f = .5), span = 1)+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))
# coord_trans(y = "log")

p.cov + p.nuc.div + plot_layout(nrow = 2, guides = "collect")






ggplot(data = one.genome, aes(x = yday, y = divergent_site_count))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
  geom_point(aes(color = season), alpha = .2)+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))

tot.selection <- one.genome[!is.na(dNdS_substitutions) & dNdS_substitutions!=0, .(.N, season,yday,year), by = .(genome, date)]

tot.selection <- one.genome[!is.na(divergent_site_count), .(.N, season,yday,year), by = .(genome, date)]

ggplot(data = tot.selection, aes(x = date, y = N))+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
  geom_point(aes(color = season), size = 3)+
  geom_line()+
  scale_x_continuous(name = element_blank(), breaks = parse_date_time(seq.int(2000,2020,5),"y"), labels = seq.int(2000,2020,5))

ggplot(data = tot.selection, aes(x = yday, y = N))+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"))+
  geom_point(aes(color = season, group = year), size = 3)+
  # geom_line(aes(group = year))+
  scale_x_continuous(name = element_blank(), breaks = yday(parse_date_time(seq.int(1:12),"m")), labels = month(parse_date_time(seq.int(1:12),"m"), abbr = T, label = T))
# no patterns there...? I don't understand what these things are anyway

tot.selection <- one.genome[!is.na(dNdS_substitutions) & dNdS_substitutions>1, .(.N, season,yday,year), by = .(genome, date)]

# OK, how would you summarize by cogs under selection?
summary <- genes[!is.na(dNdS_substitutions) & dNdS_substitutions > 1, .(.N), by = .(genome)] # and add "COG" column to by

ggplot(summary, aes(x = N, y = genome))+
  geom_bar()

