



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

