library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidystats)
library(car)

# DATA ####
setwd("CocktailvsGeneralist/Manuscript/Figure S2")

# Specify colors for each Treatment
colores <- c('#167D0A','#6038CE')
names(colores) <- c('2811','Cocktail')

# Experiment comparing 2811 and Cocktail
rawdat <- read.csv(file = "21-11-16-coevol-counts.csv")
#rawdat <- rawdat[!(rawdat$Treatment=="2811" & rawdat$Rep == '2'),]
rawdat$TxR <- paste(rawdat$Treatment, rawdat$Rep, sep = '-')
# Order Treatments for rawdat (plotting order)
rawdat$Treatment <- as.factor(rawdat$Treatment)
rawdat$Treatment <- factor(rawdat$Treatment, levels=c('2811','Cocktail'))

# SPEC INDEX ####
# Creates a plot of how specialization index changes over time
rawdat$Group <- paste(rawdat$Treatment, rawdat$Rep, rawdat$Timepoint, sep = '-')
length(unique(rawdat$Group))
datframe <- data.frame(Treatment = as.character(),
                        Rep = as.character(),
                        Timepoint = as.character(),
                        SIndex = as.character())
# limit of detection = 10000
keep_groups <- subset(rawdat, Plated == 'Phage_606')$Group
rawdat <- subset(rawdat, Group %in% keep_groups)

# Calc spec index ####
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
colnames(datframe) <- c('Timepoint','Treatment','Rep','SIndex')
datframe$Treatment <- as.factor(datframe$Treatment)
datframe$Treatment <- factor(datframe$Treatment, levels=c('2811','Cocktail'))
datframe <- subset(datframe, SIndex != 'NaN')
datframe <- datframe[order(datframe$Rep),]

# Make SI plots
lit2811 <- list()
for (i in 1:length(unique(datframe$Rep))) {
  print(i)
  R <- subset(datframe, Treatment == '2811' &
                Rep == unique(datframe$Rep)[i])
  lit2811[[i]] <- ggplot(data = R,
                         aes(x = as.numeric(Timepoint), y = SIndex,
                             color = Treatment)) +
    geom_path(size = 1) + geom_point() +
    labs(y = 'Specialization Index', x = 'Day') +
    scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1), limits = c(-1,1)) +
    scale_x_continuous(breaks = 0:16, limits = c(0,16)) +
    scale_color_manual(values = colores,
                       labels = c(expression(lambda*tgen),
                                  'Cocktail')) +
    theme_bw() +
    theme(legend.title = element_blank())
}

litC <- list()
for (i in 1:length(unique(datframe$Rep))) {
  R <- subset(datframe, Treatment == 'Cocktail' &
                Rep == unique(datframe$Rep)[i])
  litC[[i]] <- ggplot(data = R,
                         aes(x = as.numeric(Timepoint), y = SIndex,
                             color = Treatment)) +
    geom_path(size = 1) + geom_point() +
    labs(y = 'Specialization Index', x = 'Day') +
    scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1), limits = c(-1,1)) +
    scale_x_continuous(breaks = 0:16, limits = c(0,16)) +
    scale_color_manual(values = colores,
                       labels = c(expression(lambda*tgen),
                                  'Cocktail')) +
    theme_bw() +
    theme(legend.title = element_blank())
}
lit2811[[3]]
litC[[1]]
SI_matrix <- ggarrange(lit2811[[3]],litC[[1]],litC[[3]],litC[[5]],
                       lit2811[[6]],litC[[2]],litC[[4]],litC[[6]],
                         legend = 'top',
                         common.legend = T,
                         ncol=4, nrow=2)
ggsave(filename = 'SI_matrix.pdf', SI_matrix,
       width = 15, height = 7)











