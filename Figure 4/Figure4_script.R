
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidystats)
library(car)

setwd("CocktailvsGeneralist/Manuscript/Figure 4")
# Experiment comparing 2811 and Cocktail
rawdat <- read.csv(file = "21-11-16-coevol-counts.csv")
rawdat <- rawdat[!(rawdat$Treatment=="2811" & rawdat$Rep == '2'),]
rawdat$TxR <- paste(rawdat$Treatment, rawdat$Rep, sep = '-')
# Order Treatments for rawdat (plotting order)
rawdat$Treatment <- as.factor(rawdat$Treatment)
rawdat$Treatment <- factor(rawdat$Treatment, levels=c('2811','Cocktail'))

# Specify colors for each Treatment
colores <- c('#167D0A','#6038CE')
names(colores) <- c('2811','Cocktail')

# Summarize data ####
sumdat <- subset(rawdat, Plated == 'Bacteria' &
                     Titer != 'NA')  %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    mean = mean(Titer),
    sd = sd(Titer),
    n = n(),
    se = sd / sqrt(n)
  )
CFUdat <- subset(rawdat, Plated == 'Bacteria' & Titer !='NA')

# Plot ####
CFUplot <- ggplot(CFUdat, aes(x = Timepoint, y = log10(Titer),
                              color = Treatment,
                              group = interaction(Treatment, Rep),
                              linetype = TxR)) +
  geom_path(size = 0.5) +
  geom_point(size = 2) +
  labs(y = bquote('Bacteria Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
  scale_y_continuous(breaks = 0:10,
                     limits = c(-0.2,10.2),
                     expand = c(0,0)) +
  scale_x_continuous(breaks = 0:16,
                     limits = c(-0.2,16.5),
                     expand = c(0,0)) +
  scale_color_manual(values = colores,
                     labels = c(expression(lambda*tgen),
                                'Cocktail')) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(linetype = 'none') +
  geom_hline(yintercept = 1, linetype = 'dashed') +  # limit of detection
  annotate("rect",xmin=c(-0.1,-0.1),xmax=c(16.4,16.4),
           ymin=c(-0.1,-0.1) ,ymax=c(0.9,0.9),
           alpha=1, color="gray", fill="gray") +
  scale_linetype_manual(values = c(1,2,1,1,3,
                                   1,1,1,1,1,1)) +
  annotate(geom = 'text', x=1, y=9.7, label = '**') +
  annotate(geom = 'text', x=2, y=9.7, label = '**') +
  annotate(geom = 'text', x=3, y=9.7, label = '**') +
  annotate(geom = 'text', x=4, y=9.7, label = '*') +
  annotate(geom = 'text', x=5, y=9.7, label = '*') +
  annotate(geom = 'text', x=6, y=9.7, label = '*') +
  annotate(geom = 'text', x=7, y=9.7, label = '*') +
  annotate(geom = 'text', x=8, y=9.7, label = '*') +
  annotate(geom = 'text', x=9, y=9.7, label = '*') +
  annotate(geom = 'text', x=10, y=9.7, label = '*') +
  annotate(geom = 'text', x=11, y=9.7, label = '*') +
  annotate(geom = 'text', x=12, y=9.7, label = '*') +
  annotate(geom = 'text', x=13, y=9.7, label = '*') +
  annotate(geom = 'text', x=14, y=9.7, label = '*') +
  annotate(geom = 'text', x=15, y=9.7, label = '*') +
  annotate(geom = 'text', x=16, y=9.8, label = 'ns') +
  theme(axis.title = element_text(size = 12))
print(CFUplot)
ggsave('CFUplot.pdf',
       CFUplot, width = 6, height = 5)

# Statistics ####
# Initial phage titer ####
# phage.initial <- subset(rawdat, Plated == 'Phage_606' & Timepoint == '0')
# t.test(Titer~Treatment, data = phage.initial)
# 2811  lower titer (1.3*10^5) than Cocktail (8.75*10^5)

# 2811 vs Cocktail ####
vec <- vector(length = 16)
for (i in 1:16) {
  x = subset(rawdat, Timepoint == i & Plated == 'Bacteria' & Titer != 'NA')
  print(i)
  y = wilcox.test(data=x, Titer~Treatment)$p.value
  print(y)
  vec[i] <- y
}
vec
# p-values
# day 1 -- 0.004 **
# day 2 -- 0.004 **
# day 3 -- 0.004 **
# day 4 -- 0.013 *
# day 5 -- 0.008 *
# day 6 -- 0.021 *
# day 7 -- 0.008 *
# day 8 -- 0.013 *
# day 9 -- 0.013 *
# day 10 -- 0.013 *
# day 11 -- 0.011 *
# day 12 -- 0.008 *
# day 13 -- 0.008 *
# day 14 -- 0.008 *
# day 15 -- 0.021 *
# day 16 -- 0.053 ns
# significantly more 2811 suppression for 15 d

