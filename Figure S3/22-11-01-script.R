
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(car)

setwd("CocktailvsGeneralist/Manuscript/Figure S3")
rawdat <- read.csv(file = "22-11-01-data.csv")
rawdat <- subset(rawdat, Pop != '3')
rawdat <- subset(rawdat, Pop != '4')

# Manually specify colors for each Treatment
colores <- c('#66CCEE','#EE6677',
             '#F28E2B','#4E79A7',
             '#B07AA1',
             '#59A14F','#EDC948')
names(colores) <- c('Lspec','Ospec',
                    'P1T16_F-', 'P1T16_L-',
                    'P2T16_F-',
                    'P5T16_F-', 'P5T16_L-')

# CALCULATE EOP ####
# with respect to controls measured on same day
# make index and empty EOP columns
rawdat$index <- 1:nrow(rawdat)
rawdat$EOP <- 0

# subset for each phage
#  take WT titer
#  calculate EOP and append to rawdat$EOP by index
for (i in unique(rawdat$Date)) {
  for (j in unique(rawdat$Phage)) {
    # can set both vals and subset by both at same time
    sub <- subset(rawdat, Date == i & Phage == j)
    titer_WT <- sub[sub$Host == 'WT',]$Titer
    sub$EOP <- sub$Titer/titer_WT
    rawdat$EOP[sub$index] <- sub$EOP
  }
}

# log10 Trans ####
rawdat$logEOP <- log10(rawdat$EOP)

# Plot  ####
man_dat <- subset(rawdat, Host == 'ΔManY' | Host == 'ΔManZ')
ord <- c('ΔManY','ΔManZ')
man_dat$Host <- factor(man_dat$Host, levels = ord)

man_eop <- ggplot(data = man_dat, 
                  aes(x = Host, y = log10(EOP),
                      color = Phage, fill = Phage)) +
  geom_bar(stat = 'identity', position='dodge') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top') +
  labs(y = bquote(~log[10]~'(EOP)')) +
  scale_y_continuous(limits = c(-5,0.5),
                     breaks = seq(-5,0,1)) +
  scale_color_manual(values = colores,
                     guide = 'none') +
  scale_fill_manual(values = colores) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(axis.title = element_text(size = 12))
man_eop
ggsave(filename = 'mannose-eop-log.pdf', man_eop,
       width = 6, height = 6)

