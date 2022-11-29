library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidystats)
library(car)

# DATA ####
setwd("CocktailvsGeneralist/Manuscript/Figure S1")

# Specify colors for each Treatment
colores <- c('#DCBE26','#6038CE')
names(colores) <- c('EvoC','L+O')

rawdat <- read.csv(file = "2021-05-18_CoevolutionCounts.csv")
# Order Treatments for rawdat (plotting order)
rawdat = subset(rawdat, Treatment %in% c('EvoC','L+O'))
rawdat$Treatment <- as.factor(rawdat$Treatment)
rawdat$Treatment <- factor(rawdat$Treatment, levels=c('EvoC','L+O'))

# BACTERIAL TITER ####
# This set of code loops through replicate flasks for each treatment
# Then makes and saves plots into a list
loopdat <- subset(rawdat,
                  Plated == 'Bacteria')
lit <- list()
for (j in 1:length(unique(loopdat$Treatment))) {
  treats <- unique(loopdat$Treatment)[j]
  y1 <- subset(loopdat, Treatment == treats)
  if (treats == "EvoC") {
    for (w in 1:length(unique(y1$Rep))) {
      lit[[w]] <- ggplot(data = subset(y1, Rep == unique(y1$Rep)[w] & Titer !='NA'),
                         aes(x = Timepoint,
                             y = log10(as.numeric(Titer)),
                             color = Treatment)) +
        geom_path(size = 1) + geom_point() +
        labs(y = bquote('Bacterial Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
        scale_y_continuous(breaks = 0:10,limits = c(0,10),expand = c(0,0)) +
        scale_x_continuous(breaks = 0:10,limits = c(0,10),expand = c(0,0)) +
        scale_color_manual(values = colores,
                           labels = c(expression(lambda*egen),
                                      'Cocktail')) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        # hline is the limit of detection
        geom_hline(yintercept = 1, linetype = 'dashed')
    } }
  else {
    for (p in 1:length(unique(y1$Rep))) {
      lit[[p+12]] <- ggplot(data = subset(y1, Rep == unique(y1$Rep)[p] & Titer !='NA'),
                            aes(x = Timepoint,
                                y = log10(as.numeric(Titer)),
                                color = Treatment)) +
        geom_path(size = 1) + geom_point() +
        labs(y = bquote('Bacterial Titer ('~log[10]~'CFU/mL)'), x = 'Day') +
        scale_y_continuous(breaks = 0:10,limits = c(0,10),expand = c(0,0)) +
        scale_x_continuous(breaks = 0:10,limits = c(0,10),expand = c(0,0)) +
        scale_color_manual(values = colores,
                           labels = c(expression(lambda*egen),
                                      'Cocktail')) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        # hline is the limit of detection
        geom_hline(yintercept = 1, linetype = 'dashed')
    } } } 
# bact plot matrix ####
b_titer_matrix <- ggarrange(lit[[1]],lit[[2]],lit[[4]],lit[[5]],lit[[6]],lit[[7]],
                         lit[[8]],lit[[9]],lit[[11]],lit[[12]],NaN,NaN,
                         lit[[13]],lit[[14]],lit[[15]],lit[[16]],lit[[17]],lit[[18]],
                         lit[[19]],lit[[20]],lit[[22]],lit[[23]],lit[[24]],NaN,
                         legend = 'top',
                         common.legend = T,
                         ncol=6, nrow=4)
ggsave(filename = 'Bacterial-Titer-Matrix.png', b_titer_matrix,
       width = 16, height = 12)


# SPEC INDEX ####
# Creates a plot of how specialization index changes over time
rawdat$Group <- paste(rawdat$Treatment, rawdat$Rep, rawdat$Timepoint, sep = '-')
length(unique(rawdat$Group))
datframe <- data.frame(Treatment = as.character(),
                        Rep = as.character(),
                        Timepoint = as.character(),
                        SIndex = as.character())

# limit of detection = 10000
keep_groups <- subset(rawdat, Plated == 'Phage_606' & Titer > 10000)$Group
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
colnames(datframe) <- c('Treatment','Rep','Timepoint','SIndex')
datframe$Treatment <- as.factor(datframe$Treatment)
datframe$Treatment <- factor(datframe$Treatment, levels=c('EvoC','L+O'))
datframe <- subset(datframe, SIndex != 'NaN')

# Make SI plots
litty <- list()
for (t in 1:length(levels(datframe$Treatment))) {
  z <- levels(datframe$Treatment)[t]
  zz <- subset(datframe, Treatment == z)
  if (z == "EvoC") {
    for (r in 1:length(unique(zz$Rep))) {
      vv <- subset(zz, Rep == r)
      litty[[r]] <- ggplot(data = vv,
                         aes(x = as.numeric(Timepoint),
                             y = SIndex,
                             color = Treatment)) +
        geom_path(size = 1) + geom_point() +
        labs(y = 'Specialization Index', x = 'Day') +
        scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1),
                           limits = c(-1,1)) +
        scale_x_continuous(breaks = 0:10,
                           limits = c(0,10)) +
        scale_color_manual(values = colores) +
        theme_bw()
    } }
  else {
    for (y in 1:length(unique(zz$Rep))) {
      uu <- subset(zz, Rep == y)
      litty[[y+12]] <- ggplot(data = uu,
                            aes(x = as.numeric(Timepoint),
                                y = SIndex,
                                color = Treatment)) +
        geom_path(size = 1) + geom_point() +
        labs(y = 'Specialization Index', x = 'Day') +
        scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1),
                           limits = c(-1,1)) +
        scale_x_continuous(breaks = 0:10,
                           limits = c(0,10)) +
        scale_color_manual(values = colores) +
        theme_bw()
    }
  }
}
SI_matrix <- ggarrange(litty[[1]],litty[[2]],litty[[4]],litty[[5]],litty[[6]],litty[[7]], 
                         NaN,NaN,litty[[8]],litty[[9]],litty[[11]],litty[[12]],
                         litty[[13]],litty[[14]],litty[[15]],litty[[16]],litty[[17]],litty[[18]],
                         NaN,litty[[19]],litty[[20]],litty[[22]],litty[[23]],litty[[24]],
                         legend = 'top',
                         common.legend = T,
                         ncol=6, nrow=4)
ggsave(filename = 'SI_matrix.png', SI_matrix,
       width = 16, height = 12)


b_titer_matrix <- ggarrange(lit[[1]],lit[[2]],lit[[4]],lit[[5]],lit[[6]],lit[[7]],
                            lit[[8]],lit[[9]],lit[[11]],lit[[12]],NaN,NaN,
                            lit[[13]],lit[[14]],lit[[15]],lit[[16]],lit[[17]],lit[[18]],
                            lit[[19]],lit[[20]],lit[[22]],lit[[23]],lit[[24]],NaN,
                            legend = 'top',
                            common.legend = T,
                            ncol=6, nrow=4)

# combined titer and SI ####
combine_matrix <- ggarrange(lit[[1]],lit[[2]],lit[[4]],lit[[5]],lit[[6]],lit[[7]],
                            lit[[8]],lit[[9]],lit[[11]],lit[[12]],NaN,NaN,
                            litty[[1]],litty[[2]],litty[[4]],litty[[5]],litty[[6]],litty[[7]],
                            litty[[8]],litty[[9]],litty[[11]],litty[[12]],NaN,NaN,
                            lit[[13]],lit[[14]],lit[[15]],lit[[16]],lit[[17]],lit[[18]],
                            lit[[19]],lit[[20]],lit[[22]],lit[[23]],lit[[24]],NaN,
                            litty[[13]],litty[[14]],litty[[15]],litty[[16]],litty[[17]],litty[[18]],
                            litty[[19]],litty[[20]],litty[[22]],litty[[23]],litty[[24]],NaN,
                            legend = 'top',
                            common.legend = T,
                            ncol=12, nrow=4)
ggsave(filename = 'combined_matrix.png', combine_matrix,
       width = 24, height = 11)














