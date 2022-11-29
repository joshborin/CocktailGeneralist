
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(car)


setwd("CocktailvsGeneralist/Manuscript/Figure 6")
rawdat <- read.csv(file = "2022-02-12-mannose-eop.csv")

# Manually specify colors for each Treatment
colores <- c('#DCBE26','#66CCEE','#EE6677','#167D0A')
names(colores) <- c('EvoC','LamB Spec','OmpF Spec','2811')

# MANNOSE EOP assay ####
#calculate eop ####
# with respect to controls measured on same day
# make index and empty EOP columns
rawdat$index <- 1:nrow(rawdat)
rawdat$EOP <- 0

# subset for each day batch and phage
#  take REL606 titer for that day
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

# log10 transform EOP values
rawdat$logEOP <- log10(rawdat$EOP)
#rawdat$logEOP[rawdat$logEOP == -Inf] <- -3

# plot  ####
ord <- c('WT','ΔManY','ΔManZ')
rawdat$Host <- factor(rawdat$Host, levels = ord)
ord2 <- c('2811','LamB Spec','OmpF Spec','EvoC')
rawdat$Phage <- factor(rawdat$Phage, levels = ord2)

man_eop <- ggplot(data = subset(rawdat, Host != 'WT'), 
                  aes(x = Host, y = log10(EOP),
                      color = Phage, fill = Phage)) +
  geom_bar(stat = 'identity', position='dodge') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top') +
  labs(y = bquote(~log[10]~'(EOP)')) +
  scale_y_continuous(limits = c(-5,0),
                     breaks = seq(-5,0,1)) +
  scale_color_manual(values = colores,
                     guide = 'none') +
  scale_fill_manual(values = colores,
                    labels = c(expression(lambda*egen),
                               expression(lambda*Lspec),
                               expression(lambda*Ospec),
                               expression(lambda*tgen)),
                    guide = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(axis.title = element_text(size = 12))
man_eop
ggsave(filename = 'mannose-eop-log-2.pdf', man_eop,
       width = 4, height = 5)














