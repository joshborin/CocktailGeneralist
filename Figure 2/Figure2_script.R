
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidystats)
library(car)

# Specify colors for each Treatment
colores <- c('black','#DCBE26','#6038CE')
names(colores) <- c('cI26','EvoC','L+O')

# First experiment with more treatment comparisons, fewer replicates
rawdat2 <- read.csv(file = "2021-03-28_CoevolutionCounts.csv")
# Order Treatments for rawdat2 (plotting order)
rawdat2 <- subset(rawdat2, Treatment == 'cI26' |
                    Treatment == 'EvoC' |
                    Treatment == 'L+O')
rawdat2$Treatment <- as.factor(rawdat2$Treatment)
rawdat2$Treatment <- factor(rawdat2$Treatment, levels=c('cI26','EvoC','L+O'))

# Second, larger experiment comparing EvoC and Cocktail
rawdat <- read.csv(file = "2021-05-18_CoevolutionCounts.csv")
# Order Treatments for rawdat (plotting order)
rawdat$Treatment <- as.factor(rawdat$Treatment)
rawdat$Treatment <- factor(rawdat$Treatment, levels=c('EvoC','L+O'))

# EXPERIMENT 1 ####
# PANEL A ####
# sum data ####
sumdat2 <- subset(rawdat2, Plated == 'Bacteria' &
                     Titer != 'NA')  %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    mean = mean(Titer),
    sd = sd(Titer),
    n = n(),
    se = sd / sqrt(n),
    medi = median(Titer)
  )
CFUdat2 <- subset(rawdat2, Plated == 'Bacteria' & Titer !='NA')

# plot data ####
PanelA <- ggplot(CFUdat2, 
                 aes(x = Timepoint, y = log10(Titer),
                     color = Treatment,
                     group = interaction(Treatment, Rep),
                     shape = as.factor(Rep))) +
  geom_line(size = 0.5) +
  geom_point(size = 2) +
  labs(y = bquote('Bacteria Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = 0:10,
                     limits = c(-0.2,10.2),
                     expand = c(0,0)) +
  scale_x_continuous(breaks = 0:10,
                     limits = c(-0.2,10.5),
                     expand = c(0,0)) +
  scale_color_manual(values = colores,
                     labels = c(expression(lambda*anc),
                                expression(lambda*egen),
                                'Cocktail')) +
  scale_shape_discrete(guide = 'none') +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  geom_hline(yintercept = 1, linetype = 'dashed') + # limit of detection
  annotate("rect", xmin=c(-0.1,-0.1), xmax=c(10.4,10.4),
           ymin=c(-0.1,-0.1) , ymax=c(0.9,0.9), alpha=1, color="gray", fill="gray") +
  theme(axis.title = element_text(size = 12))
print(PanelA)
ggsave('PanelA.png',
       PanelA, width = 6, height = 4)

# statistics ####
# Initial phage titer ####
#phage.initial <- subset(rawdat2, Plated == 'Phage_606' & Timepoint == '0' &
#                          Treatment != 'cI26')
#t.test(Titer~Treatment, data = phage.initial)
# L+O cocktail lower titer (1.1*10^6) than EvoC (1.6*10^6)

# cI26 vs EvoC ####
vec <- vector(length = 10)
for (i in 1:10) {
  x = subset(rawdat2, Timepoint == i & Plated == 'Bacteria' & Titer != 'NA' &
               Treatment != 'L+O')
  print(i)
  y = wilcox.test(data=x, Titer~Treatment, alternative = 't')$p.value
  print(y)
  vec[i] <- y
}
vec
# p-values
# day 1 -- 0.1 
# day 2 -- 0.1 
# day 3 -- 0.7 
# day 4 -- 0.2 
# day 5 -- 0.1 
# day 6 -- 0.1 
# day 7 -- 0.08
# day 8 -- 0.1 
# day 9 -- 0.2 
# day 10 -- 0.4 
# significant for first 2 days at p=0.1
# note that if computing a one-tailed test, the minimum p-val is 0.05
# but with limited n =3, the minimum for a two-tailed is p=0.1

# cI26 vs L+O ####
vec <- vector(length = 10)
for (i in 1:10) {
  x = subset(rawdat2, Timepoint == i & Plated == 'Bacteria' & Titer != 'NA' &
               Treatment != 'EvoC')
  print(i)
  y = wilcox.test(data=x, Titer~Treatment)$p.value
  print(y)
  vec[i] <- y
}
vec
# p-values
# day 1 -- 0.028 
# day 2 -- 0.024 
# day 3 -- 0.154 
# day 4 -- 0.143 
# day 5 -- 0.048 
# day 6 -- 0.154 
# day 7 -- 0.092
# day 8 -- 0.052 
# day 9 -- 0.897 
# day 10 -- 0.517 
# significant for first 2 days at p=0.05

# EvoC vs L+O ####
vec <- vector(length = 10)
for (i in 1:10) {
  x = subset(rawdat2, Timepoint == i & Plated == 'Bacteria' & Titer != 'NA' &
               Treatment != 'cI26')
  print(i)
  y = wilcox.test(data=x, Titer~Treatment)$p.value
  print(y)
  vec[i] <- y
}
vec
# p-values
# day 1 -- 0.154 
# day 2 -- 0.048 
# day 3 -- 0.243 
# day 4 -- 0.857 
# day 5 -- 0.905 
# day 6 -- 1.00 
# day 7 -- 0.345
# day 8 -- 0.694 
# day 9 -- 0.511 
# day 10 -- 0.896 
# significant one day 2, otherwise the same


# EXPERIMENT 2 ####
# omit flasks where phage stochastically lost (cells not resistant)
evoc_dat <- subset(rawdat, Treatment %in% 'EvoC' & Rep != '3'
                   & Rep != '10')
LO_dat <- subset(rawdat, Treatment %in% 'L+O' & Rep != '9')
rawdat <- rbind(evoc_dat, LO_dat)

# PANEL B ####
# sum data ####
sumdat <- subset(rawdat, Plated == 'Bacteria' &
                    Titer != 'NA')  %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    mean = mean(Titer),
    sd = sd(Titer),
    n = n(),
    se = sd / sqrt(n),
    medi = median(Titer))
CFUdat <- subset(rawdat, Plated == 'Bacteria' & Titer !='NA')

#plot data ####
PanelB <- ggplot(sumdat, aes(x = Timepoint, y = log10(medi), color = Treatment)) +
  geom_line(data = CFUdat, aes(x = Timepoint, y = log10(Titer),
                               group = interaction(Treatment, Rep)),
            size = 0.8, alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  labs(y = bquote('Bacteria Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = 0:10,
                     limits = c(-0.2,10.2),
                     expand = c(0,0)) +
  scale_x_continuous(breaks = 0:10,
                     limits = c(-0.2,10.5),
                     expand = c(0,0)) +
  scale_color_manual(values = colores,
                     labels = c(expression(lambda*anc),
                                expression(lambda*egen),
                                'Cocktail')) +
  theme_bw(base_size = 16) +
  geom_hline(yintercept = 1, linetype = 'dashed') + # limit of detection
  annotate("rect", xmin=c(-0.1,-0.1), xmax=c(10.4,10.4),
           ymin=c(-0.1,-0.1) , ymax=c(0.9,0.9), alpha=1, color="gray", fill="gray") +
  annotate(geom = 'text', x=1, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=2, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=3, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=4, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=5, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=6, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=7, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=8, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=9, y=9.7, label = 'ns') +
  annotate(geom = 'text', x=10, y=9.7, label = 'ns') +
  theme(axis.title = element_text(size = 12))
print(PanelB)
ggsave('PanelB.png',
       PanelB, width = 6, height = 4)

PanelAB <- ggarrange(PanelA, PanelB,
                         legend = 'top', common.legend = TRUE,
                         ncol = 1, nrow = 2, labels = c('A','B'))
ggsave(filename = 'PanelsAB-2.pdf', PanelAB,
       width = 5, height = 7)

# statistics ####
# Initial phage titers ####
#phage.initial <- subset(rawdat, Plated == 'Phage_606' & Timepoint == '0')
#t.test(Titer~Treatment, data = phage.initial)
# L+O cocktail higher titer (5.5*10^5) than EvoC (1.2*10^5)

# Diff suppress ####
vec <- vector(length = 10)
for (i in 1:10) {
  x = subset(rawdat, Timepoint == i & Plated == 'Bacteria' & Titer != 'NA')
  print(i)
  y = wilcox.test(data=x, Titer~Treatment)$p.value
  print(y)
  vec[i] <- y
}
vec
# p-values
# day 1 -- 0.398 
# day 2 -- 0.809 
# day 3 -- 0.113 
# day 4 -- 0.863 
# day 5 -- 0.756 
# day 6 -- 1.000 
# day 7 -- 0.705
# day 8 -- 0.349 
# day 9 -- 0.673 
# day 10 -- 0.223 

