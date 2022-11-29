
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(car)

setwd("CocktailvsGeneralist/2022-11-01-evo-cocktail")
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

# SI INDEX ####
#rawdat$Group <- paste(rawdat$Phage, rawdat$Rep, rawdat$Timepoint, sep = '_')
length(unique(rawdat$Phage))
datframe <- data.frame(Phage = as.character(),
                       SIndex = as.character())

for (i in 1:length(unique(rawdat$Phage))) {
  print(i)
  print(unique(rawdat$Phage)[i])
  x <- subset(rawdat, Phage == unique(rawdat$Phage)[i])
  L <- subset(x, Host == 'ΔOmpF')$Titer
  O <- subset(x, Host == 'ΔLamB')$Titer
  dat <- data.frame(nrow = 1, ncol = 2)
  dat[1,1] <- as.character(x[1,5])
  dat[1,2] <- ((L-O)/(L+O))
  datframe <- rbind(as.data.frame(dat), datframe)
}
colnames(datframe) <- c('Phage','SIndex')
datframe$Phage <- as.factor(datframe$Phage)
levels(datframe$Phage)

# Plot ####
SIplot <- ggplot(data = datframe, 
                  aes(x = Phage, y = SIndex,
                      color = Phage, fill = Phage)) +
  geom_bar(stat = 'identity', position='dodge') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top') +
  labs(y = bquote(~log[10]~'(EOP)')) +
  scale_y_continuous(limits = c(-1,1),
                     breaks = seq(-1,1,0.5)) +
    scale_color_manual(values = colores,
                       guide = 'none') +
    scale_fill_manual(values = colores) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(axis.title = element_text(size = 12)) +
  coord_flip()
SIplot
ggsave(filename = 'SIndex.pdf', SIplot,
       width = 8, height = 4)




















