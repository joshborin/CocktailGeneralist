library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(car)

setwd("CocktailvsGeneralist/Manuscript/Figure 5")
rawdat <- read.csv(file = "22-01-20-resistance-master.csv")
rawdat <- rawdat[!is.na(rawdat$Titer),]

# Manually specify colors for each Treatment
colores <- c('#DCBE26','#66CCEE','#EE6677','#167D0A')
names(colores) <- c('EvoC','LamB Spec','OmpF Spec','2811')

# CALCULATE EOP ####
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
    titer_606 <- sub[sub$Treatment == 'REL606',]$Titer
    sub$EOP <- sub$Titer/titer_606
    rawdat$EOP[sub$index] <- sub$EOP
  }
}

# log10 Trans ####
rawdat$logEOP <- log10(rawdat$EOP)
#rawdat$logEOP[rawdat$logEOP == -Inf] <- -8

# Log10 line plot EOP ####
sumdat <- subset(rawdat, Population != 0) %>%
  group_by(Date,Treatment,Population,Timepoint,Phage) %>%
  summarise(
    mean = (mean(EOP)),
    median = (median(EOP)),
    sd = (sd(EOP)),
    n = n(),
    se = sd / sqrt(n) 
  )

sumdat$logmedi <- log10(sumdat$median)
#sumdat$logmedi[sumdat$logmedi == -Inf] <- -8

pd <- position_dodge(0.3)
eop_line_plot <- list()
counter2 <- 0
for (i in 1:length(unique(sumdat$Treatment))) {
  treat <- unique(sumdat$Treatment)[i]
  sum_T <- subset(sumdat, Treatment == treat)
  sub_T <- subset(rawdat, Treatment == treat)
  for (j in 1:length(unique(sort(sub_T$Population)))) {
    pop <- unique(sort(sub_T$Population))[j]
    if (treat == 'REL606' | treat == "L-" | treat == 'O-' | pop < 1) {
    } else {
      counter2 <- counter2 + 1
      sum_TP <- subset(sum_T, Treatment == treat &
                         Population == pop)
      sub_TP <- subset(sub_T, Treatment == treat &
                         Population == pop)
      eop_line_plot[[counter2]] <- ggplot(data = sum_TP, 
                                          aes(x = Timepoint,
                                              y = log10(median),
                                              color = Phage)) +
        theme_classic() +
        geom_line(position = pd, size = 1.5, alpha = 1) +
        geom_point(position = pd, size = 2, alpha = 1) +
        #geom_point(data = sub_TP,
        #            aes(x = Timepoint, y = log10(EOP)),
        #            position = pd, alpha = 0.1) +
        labs(title = paste(treat,'Population',pop,sep = ' '),
             x = 'Day',
             y = bquote(~log[10]~'(EOP)')) +
        scale_y_continuous(limits = c(-8.3,1),
                           breaks = seq(-8,0,1),
                           expand = c(0,0)) +
        scale_x_continuous(limits = c(-0.5,16.5),
                           breaks = c(0,1,2,3,4,5,10,16),
                           expand = c(0,0)) +
        scale_color_manual(values = colores,
                           labels = c(expression(lambda*egen),
                                      expression(lambda*Lspec),
                                      expression(lambda*Ospec),
                                      expression(lambda*tgen))) +
        scale_fill_manual(values = colores) +
        geom_hline(yintercept = 0, linetype = 'dashed')
        #annotate("rect",xmin=c(-0.5,-0.5),xmax=c(16.5,16.5),
        #         ymin=c(-8.3,-8.3) ,ymax=c(-8,-8),
        #         alpha=1, color="gray", fill="gray")
    }
  }
}
eop_line_plot[[4]]
EOP_plot_mat_tog <- ggarrange(eop_line_plot[[8]], eop_line_plot[[1]], eop_line_plot[[2]],eop_line_plot[[3]],
                              eop_line_plot[[9]], eop_line_plot[[4]], eop_line_plot[[5]],eop_line_plot[[6]],
                              legend = 'top',
                              common.legend = TRUE,
                              ncol =4, nrow =2)
#EOP_plot_mat_tog
ggsave(filename = 'logEOP-median2.pdf', EOP_plot_mat_tog,
       width = 12, height = 5)


