# RRR
# same old plot from PNAS. getting a lot of use... hope that's ok :|
# fig 1a

library(data.table)
library(lubridate)
# source(file = "scripts-figures/2024-01-01_so_much_data_figure/plot_sample_dates_fancy.R")
source(file = "pop/scripts-generic_plots/plot_sample_dates_fancy.R")

seas <- readRDS(file = "figures/2023-12-24_paper_figures/Rohwer_Figure_1_data/season_dates.rds")

sample.key <- fread(file = "data/2023-12-08_combined_genome_info/sample_key.tsv", colClasses = c("date" = "character"))
sample.key[ ,date := parse_date_time(x = date, orders = "ymd", tz = "Etc/GMT-5")]

# season shading 

set.up.plot(dates.vector = sample.key$date, year.lines = F, month.lines = F, 
            x.ax.lab.line = -.65, x.ax.single.letters = F,
            y.ax.lab.line = -.35, y.ax.btwn.lines = F, x.ax.btwn.lines = F,
            y.ax.tick.freq = 3, x.ax.tick.freq = 1)
# mtext(text = "Sample Year", side = 2, line = 1.9, cex = 1)
# mtext(text = "Sample Month", side = 1, line = .75, cex = 1)


shade.by.season(season.starts = seas$Spring, season.ends = seas$Clearwater, season.color = adjustcolor("tan4",.6))

shade.by.season(season.starts = seas$Clearwater, season.ends = seas$`Early Summer`, season.color = adjustcolor("cornflowerblue",.6))

shade.by.season(season.starts = seas$`Early Summer`, season.ends = seas$`Late Summer`, season.color =  adjustcolor("chartreuse4",.6))

shade.by.season(season.starts = seas$`Late Summer`, season.ends = seas$Fall, season.color = adjustcolor("purple",.4))

shade.by.season(season.starts = seas$Fall, season.ends = seas$fall.fall.end, season.color = adjustcolor("hotpink2",.5))
shade.by.season(season.starts = seas$spring.fall.start, season.ends = seas$spring.ice.on, season.color = adjustcolor("hotpink2",.5))

shade.by.season(season.starts = seas$fall.ice.on, season.ends = seas$fall.ice.off, season.color = adjustcolor("snow3",.3))
shade.by.season(season.starts = seas$spring.ice.on, season.ends = seas$Spring, season.color =  adjustcolor("snow3",.3))


add.dates(date.vector = sample.key$date, year.labs = 2000:2019, 
          point.type = "|", point.col = adjustcolor("black", alpha.f = .2), line.col = adjustcolor("black",1), 
          point.cex = .7)
