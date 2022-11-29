library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidystats)
library(car)

# DATA ####
setwd("CocktailvsGeneralist/Manuscript/Figure 3")

# Specify colors for each Treatment
colores <- c('#DCBE26','#6038CE')
names(colores) <- c('EvoC','L+O')

# EXPERIMENT 1 ####
rawdat <- read.csv(file = "2021-03-28_CoevolutionCounts.csv")
# Order Treatments for rawdat (plotting order)
rawdat = subset(rawdat, Treatment %in% c('EvoC','L+O'))
rawdat$Treatment <- as.factor(rawdat$Treatment)
rawdat$Treatment <- factor(rawdat$Treatment, levels=c('EvoC','L+O'))

# SPEC INDEX ####
# Creates a plot of how specialization index changes over time
rawdat$Group <- paste(rawdat$Treatment, rawdat$Rep, rawdat$Timepoint, sep = '_')
length(unique(rawdat$Group))
datframe <- data.frame(Treatment = as.character(),
                       Rep = as.character(),
                       Timepoint = as.character(),
                       SIndex = as.character())
# Calculate specialization index
for (i in 1:length(unique(rawdat$Group))) {
  print(i)
  print(unique(rawdat$Group)[i])
  x <- subset(rawdat, Group == unique(rawdat$Group)[i])
  L <- subset(x, Plated == 'Phage_Ominus')$Titer
  O <- subset(x, Plated == 'Phage_Lminus')$Titer
  dat <- data.frame(nrow = 1, ncol = 4)
  dat[1,1] <- as.character(x[1,2])
  dat[1,2] <- as.character(x[1,3])
  dat[1,3] <- as.character(x[1,4])
  dat[1,4] <- ((L-O)/(L+O))
  datframe <- rbind(as.data.frame(dat), datframe)
}
colnames(datframe) <- c('Treatment','Rep','Timepoint','SIndex')
datframe$Treatment <- as.factor(datframe$Treatment)
datframe$Treatment <- factor(datframe$Treatment, levels=c('EvoC','L+O'))

# SI PLOTS ####
# EvoC ####
evoCplot <- ggplot(data = subset(datframe, Treatment == 'EvoC'),
       aes(x = as.numeric(Timepoint),
           y = SIndex,
           group = Rep,
           color = Treatment,
           shape = as.factor(Rep))) +
  geom_path(size = 0.5) +
  geom_point(size = 2) +
  labs(y = 'Specialization Index', x = 'Day') +
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1),
                     limits = c(-1,1)) +
  scale_x_continuous(breaks = seq(0,10),
                     limits = c(0,10)) +
  scale_color_manual(values = colores,
                     labels = c(expression(lambda*egen),
                                'Cocktail')) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(shape = "none") +
  theme(axis.title = element_text(size = 12))
evoCplot

# Cocktail ####
cocktailplot <- ggplot(data = subset(datframe, Treatment == 'L+O'),
                   aes(x = as.numeric(Timepoint),
                       y = SIndex,
                       group = Rep,
                       color = Treatment,
                       shape = as.factor(Rep))) +
  geom_path(size = 0.5) +
  geom_point(size = 2) +
  labs(y = 'Specialization Index', x = 'Day') +
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1),
                     limits = c(-1,1)) +
  scale_x_continuous(breaks = seq(0,10),
                     limits = c(0,10)) +
  scale_color_manual(values = colores,
                     labels = c(expression(lambda*egen),
                                'Cocktail')) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(shape = 'none') +
  theme(axis.title = element_text(size = 12))
cocktailplot

# Plot Matrix ####
SIplot_mat <- ggarrange(cocktailplot,evoCplot,
                        labels = c('A','B'),
                        legend = 'top', common.legend = T,
                        ncol=1, nrow=2)
SIplot_mat
ggsave(filename = 'SpecIndex_matrix1.pdf', SIplot_mat,
       width = 5, height = 7)

