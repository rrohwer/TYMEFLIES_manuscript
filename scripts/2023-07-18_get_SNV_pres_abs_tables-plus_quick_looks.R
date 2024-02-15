# RRR
# Get SNV pres-abs tables
# goal, track when new SNVs arrive/become fixed
# but in this script just make the tables

# ---- set-up ----

library(data.table)
library(lubridate)

# local path testing
per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-21_3300043464_group3_bin69_SNVs.tsv.gz" # B
per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2016-07-20_3300033996_group7_bin32_SNVs.tsv.gz" # C
per.genome.snv.file <- "data/2023-07-16_SNVs_example_files/generated_per-genome/ME2011-09-04_3300044729_group3_bin142_SNVs.tsv.gz" # A

num.threads <- 1

output.pres.abs.folder <- "data/2023-07-16_SNVs_example_files/generated_pres_abs"

# ---- parse into pres/abs table for SNVs ----

x <- "character"
names(x) <- "date"
snv.instrain <- fread(file = per.genome.snv.file, nThread = num.threads, colClasses = x)
my.genome <- snv.instrain[1, genome]
snv.instrain[ ,is.pres := TRUE]
# snv <- snv.instrain[ ,.(scaffold, position, date)]
# snv[ ,is.pres := TRUE]

snv <- dcast(data = snv.instrain, formula = scaffold + position + gene + mutation_type ~ date, fun.aggregate = any, value.var = "is.pres", fill = FALSE)

snv[ ,snv.num := paste0("snv",1:nrow(snv))]
snv.key <- snv[ ,.(scaffold, position, snv.num, gene, mutation_type)]
snv[ ,`:=`(scaffold = NULL, position = NULL, gene = NULL, mutation_type = NULL)]
setcolorder(x = snv, neworder = c(ncol(snv), 2:(ncol(snv) - 1)))

# fwrite(x = snv, file = file.path(output.pres.abs.folder, paste0(my.genome, "_SNV_pres_abs.tsv.gz")), compress = "gzip", sep = "\t")
# fwrite(x = snv.key, file = file.path(output.pres.abs.folder, paste0(my.genome, "_SNV_key.tsv.gz")), compress = "gzip", sep = "\t")

# ---- parse into new/recurring matrix ----

# "new" SNV appearances are when it steps from 0 to 1
# subtracting offset matrices can ID this

snv.mat <- as.matrix(x = snv, rownames = "snv.num")

snv.mat <- snv.mat[ ,-1] - snv.mat[ ,-ncol(snv.mat)]

snv.mat <- snv.mat == 1

# "recurring" SNVs are when it's new but the cumulative sum of times it was new is greater than 1
snv.mat.2 <- apply(X = snv.mat, MARGIN = 1, FUN = cumsum)
snv.mat.2 <- t(snv.mat.2)
snv.mat[1:10,1:15]
snv.mat.2[1:10,1:15]

my.mat <- matrix(data = NA, nrow = nrow(snv.mat), ncol = ncol(snv.mat), dimnames = list(row.names(snv.mat), colnames(snv.mat)))
my.mat[snv.mat == TRUE & snv.mat.2 == 1] <- "N"
my.mat[snv.mat == TRUE & snv.mat.2 != 1] <- "R"

my.mat[1:10,1:15]

# ---- get number new/ number recurring per day ----

my.nonsyn <- my.mat[snv.key$mutation_type == "N", ]
my.syn  <- my.mat[snv.key$mutation_type == "S", ]

snv.stats <- data.table("Date" = colnames(my.mat), 
                        "New" = colSums(my.mat == "N", na.rm = T), 
                        "Recurring" = colSums(my.mat == "R", na.rm = T))
snv.stats[ ,`:=`(Either = New + Recurring,
                 New.Nonsyn = colSums(my.nonsyn == "N", na.rm = T),
                 Recurring.Nonsyn  = colSums(my.nonsyn == "R", na.rm = T))]
snv.stats[ ,`:=`(Either.Nonsyn = New.Nonsyn + Recurring.Nonsyn,
                 New.Syn = colSums(my.syn == "N", na.rm = T),
                 Recurring.Syn  = colSums(my.syn == "R", na.rm = T))]
snv.stats[ ,`:=`(Either.Syn = New.Syn + Recurring.Syn)]


dates.key <- unique(snv.instrain[ ,.(date, year, yday, season, invasion)])

snv.stats <- merge(x = snv.stats, y = dates.key, by.x = "Date", by.y = "date", all.x = T, all.y = F)

# ---- prep to make plots (manually for now) ----
snv.stats[ ,`:=`(Date = parse_date_time(Date, orders = "ymd"),
                 season = factor(season, levels = c("Ice.On","Spring","Clearwater","Early.Summer","Late.Summer","Fall")),
                 invasion = factor(invasion, levels = c("none","spiny","zebra")))]

acI.B <- snv.stats
acI.C <- snv.stats
acI.A <- snv.stats

acI.A <- readRDS("figures/2023-07-18_tracking_new_SNVs/acI_A.rds")
acI.B <- readRDS("figures/2023-07-18_tracking_new_SNVs/acI_B.rds")
acI.C <- readRDS("figures/2023-07-18_tracking_new_SNVs/acI_C.rds")

library(ggplot2)
library(patchwork)

# acI-B all ----

p.b.new <- ggplot(data = acI.B, aes(x = Date, y = New))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line() + 
  geom_point(aes(color = season), size = 3, alpha = .5) +
  labs(title = "acI-B", subtitle = "ME2011-09-21_3300043464_group3_bin69") +
  ylab(label = "Number new SNVs not previously observed")
p.b.new

p.b.new.trunc <- ggplot(data = acI.B, aes(x = Date, y = New))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line() + 
  geom_point(aes(color = season), size = 3, alpha = .5) +
  coord_cartesian(ylim = c(0,1000), expand = FALSE)+
  labs(title = "acI-B", subtitle = "ME2011-09-21_3300043464_group3_bin69") +
  ylab(label = "Number new SNVs not previously observed")
p.b.new.trunc

p.b.either <- ggplot(data = acI.B, aes(x = Date, y = Either))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line() + 
  geom_point(aes(color = season), size = 3, alpha = .5) + 
  labs(title = "acI-B", subtitle = "ME2011-09-21_3300043464_group3_bin69") +
  ylab(label = "Number new SNVs (including previously observed)")
p.b.either

# acI-B nonsynonymous ----

p.b.new.N <- ggplot(data = acI.B, aes(x = Date, y = New.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line() + 
  geom_point(aes(color = season), size = 3, alpha = .5) +
  labs(title = "acI-B", subtitle = "ME2011-09-21_3300043464_group3_bin69") +
  ylab(label = "Number new nonsynonymous SNVs not previously observed")
p.b.new.N

p.b.new.trunc.N <- ggplot(data = acI.B, aes(x = Date, y = New.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line() + 
  geom_point(aes(color = season), size = 3, alpha = .5) +
  coord_cartesian(ylim = c(0,500), expand = FALSE)+
  labs(title = "acI-B", subtitle = "ME2011-09-21_3300043464_group3_bin69") +
  ylab(label = "Number new nonsynonymous SNVs not previously observed")
p.b.new.trunc.N

p.b.either.N <- ggplot(data = acI.B, aes(x = Date, y = Either.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line() + 
  geom_point(aes(color = season), size = 3, alpha = .5) + 
  labs(title = "acI-B", subtitle = "ME2011-09-21_3300043464_group3_bin69") +
  ylab(label = "Number new nonsynonymous SNVs (including previously observed)")
p.b.either.N

# acI-A all ----

p.a.new <- ggplot(data = acI.A, aes(x = yday, y = New))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .1) + 
  geom_point(aes(color = season), size = 5, alpha = .5) +
  coord_cartesian(ylim = c(0,1000), expand = TRUE) + 
  labs(title = "acI-A", subtitle = "ME2011-09-04_3300044729_group3_bin142") +
  ylab(label = "Number new SNVs (previously unobserved)")
p.a.new

p.a.either <- ggplot(data = acI.A, aes(x = yday, y = Either))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .1) + 
  geom_point(aes(color = season), size = 5, alpha = .5) +
  labs(title = "acI-A", subtitle = "ME2011-09-04_3300044729_group3_bin142") +
  ylab(label = "Number new SNVs (including previously observed)")
p.a.either

p.a.new.LT <- ggplot(data = acI.A, aes(x = Date, y = New.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line() + 
  geom_point(aes(color = season), size = 3, alpha = .5) +
  labs(title = "acI-A", subtitle = "ME2011-09-21_3300043464_group3_bin69") +
  ylab(label = "Number new nonsynonymous SNVs not previously observed")
p.a.new.LT

# acI-A nonsynonymous ----

p.a.new.N <- ggplot(data = acI.A, aes(x = yday, y = New.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .1) + 
  geom_point(aes(color = season), size = 5, alpha = .5) +
  coord_cartesian(ylim = c(0,500), expand = TRUE) + 
  labs(title = "acI-A", subtitle = "ME2011-09-04_3300044729_group3_bin142") +
  ylab(label = "Number new nonsynonymous SNVs (previously unobserved)")
p.a.new.N

p.a.either.N <- ggplot(data = acI.A, aes(x = yday, y = Either.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .1) + 
  geom_point(aes(color = season), size = 5, alpha = .5) +
  labs(title = "acI-A", subtitle = "ME2011-09-04_3300044729_group3_bin142") +
  ylab(label = "Number new nonsynonymous SNVs (including previously observed)")
p.a.either.N

# acI-C all ----

p.c.new <- ggplot(data = acI.C, aes(x = yday, y = New))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .1) + 
  geom_point(aes(color = season), size = 5, alpha = .5) +
  coord_cartesian(ylim = c(0,3000), expand = TRUE) + 
  labs(title = "acI-C", subtitle = "ME2016-07-20_3300033996_group7_bin32") +
  ylab(label = "Number new SNVs (previously unobserved)")
p.c.new

p.c.either <- ggplot(data = acI.C, aes(x = yday, y = Either))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .1) + 
  geom_point(aes(color = season), size = 5, alpha = .5) +
  labs(title = "acI-C", subtitle = "ME2016-07-20_3300033996_group7_bin32") +
  ylab(label = "Number new SNVs (including previously observed)")
p.c.either

p.c.new.LT <- ggplot(data = acI.C, aes(x = Date, y = New.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line() + 
  geom_point(aes(color = season), size = 3, alpha = .5) +
  labs(title = "acI-C", subtitle = "ME2011-09-21_3300043464_group3_bin69") +
  ylab(label = "Number new nonsynonymous SNVs not previously observed")
p.c.new.LT

# acI-C nonsynonymous ----

p.c.new.N <- ggplot(data = acI.C, aes(x = yday, y = New.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .1) + 
  geom_point(aes(color = season), size = 5, alpha = .5) +
  coord_cartesian(ylim = c(0,1000), expand = TRUE) + 
  labs(title = "acI-C", subtitle = "ME2016-07-20_3300033996_group7_bin32") +
  ylab(label = "Number new nonsynonymous SNVs (previously unobserved)")
p.c.new.N

p.c.either.N <- ggplot(data = acI.C, aes(x = yday, y = Either.Nonsyn))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("snow3","tan4","cornflowerblue","chartreuse4","purple","hotpink2"), labels = c("Ice-on","Spring","Clearwater","Early Summer","Late Summer","Fall"))+
  geom_line(aes(group = year), alpha = .1) + 
  geom_point(aes(color = season), size = 5, alpha = .5) +
  labs(title = "acI-C", subtitle = "ME2016-07-20_3300033996_group7_bin32") +
  ylab(label = "Number new nonsynonymous SNVs (including previously observed)")
p.c.either.N

# combine some plots ----

p.b.new / p.a.new / p.c.new + plot_layout(guides = "collect") &
  ylab(label = "New SNVs") 

p.b.new.N / p.a.new.N / p.c.new.N + plot_layout(guides = "collect") &
  ylab(label = "New nonsynonymous SNVs") 

p.b.either / p.a.either / p.c.either + plot_layout(guides = "collect") &
  ylab(label = "New SNVs (incl. seen-before)") 

p.b.either.N / p.a.either.N / p.c.either.N + plot_layout(guides = "collect") &
  ylab(label = "New nonsynonymous SNVs (incl.seen-before") 

p.b.new / p.a.new.LT / p.c.new.LT + plot_layout(guides = "collect") &
  ylab(label = "New SNVs") 
