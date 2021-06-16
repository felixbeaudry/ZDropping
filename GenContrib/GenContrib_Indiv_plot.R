#15 June 2021
#Rose Driscoll and Felix Beaudry
#script to plot the genetic & genealogical contributions of male and female breeders for autosomes and Z
#simplified version with just the plotting commands

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(kinship2)
library(cowplot)

# load Rdata file 
setwd('~/Google Drive/Research/Data2/fsj/ZDropping_all/')

load("plotIndivGenContrib_tidy_20190529.Rdata") # code to produce all of these tables from raw data is in plotIndivGenContrib_tidy_20190425.R

#plot theme
plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3),
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"),
                    plot.background = element_rect(fill = "white"),
                    axis.text.x = element_text(size=6), axis.text.y = element_text(size=6),
                    axis.title = element_text(size=7), plot.title = element_text(size=8),
                    legend.position="right", legend.text = element_text(size=7),
                    legend.title = element_text(size=8), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(1,"cm"))


####Pop Level distributions####
#genetic (A & Z) vs. genealogical contributions in 2013 for all 926 breeders

#count of individuals with zero descendants in 2013
length(is.na(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1]))
length(is.na(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2]))

## genealogical contributions
#range and mean genealogical contributions for males with >0 descendants in 2013
min(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1]))
max(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1]))
mean(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1]))#males

#range and mean genealogical contributions for females with >0 descendants in 2013
min(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2]))
max(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2]))
mean(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2]))#females

###expected genetic contributions
##Autosomal
#range and mean expected contributions for males with >0 descendants in 2013
min(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==1 & !is.na(indiv_contribs_prop_nestlings_2013$num_descendants_2013)]))
max(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==1]))
mean(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==1]))#males

#range and mean expected contributions for females with >0 descendants in 2013
min(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==2 & !is.na(indiv_contribs_prop_nestlings_2013$num_descendants_2013)]))
max(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==2]))
mean(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==2]))#females

##Z
#males
min(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==1 & !is.na(indiv_contribs_prop_nestlings_2013$num_descendants_2013) ]))
max(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==1]))
mean(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==1]))#males

#females
min(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==2 & !is.na(indiv_contribs_prop_nestlings_2013$num_descendants_2013)]))
max(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==2]))
mean(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==2]))#females

####descents vs genes####
#genealogical vs autosomal expected model
A_FM_mod_lm <- lm(auto_mean ~  prop_2013_nestlings*Sex  , data=indiv_contribs_prop_nestlings_2013)  
summary(A_FM_mod_lm)

#genealogical vs Z expected model
Z_FM_mod_lm <- lm(Z_mean ~  prop_2013_nestlings*Sex  , data=indiv_contribs_prop_nestlings_2013)  
summary(Z_FM_mod_lm)

-0.0327305/0.0539405 #sex effect


#AZ

AZ_mod <- lm(Z_mean ~  auto_mean  , data=indiv_contribs_A_Z_2013_sex)  
summary(AZ_mod)  

AZ_mods_sex <- lm(Z_mean ~  auto_mean*Sex  , data=indiv_contribs_A_Z_2013_sex)  
summary(AZ_mods_sex)  

-5.690e-01/1.239e+00
1.239e+00 - 5.690e-01

#plot
A_contribs_vs_descendants <- 
  ggplot(indiv_contribs_prop_nestlings_2013, aes(x = prop_2013_nestlings, y = auto_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  geom_abline(slope = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Genealogical contribution in 2013", y = "Autosomal expected \ngenetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) + ylim(0,0.04)


Z_contribs_vs_descendants <- 
  ggplot(indiv_contribs_prop_nestlings_2013, aes(x = prop_2013_nestlings, y = Z_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Genealogical contribution in 2013", y = "Z expected genetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) + ylim(0,0.04)

A_vs_Z_contribs <- 
  ggplot(indiv_contribs_A_Z_2013_sex, aes(x = auto_mean, y = Z_mean)) +
    geom_abline(slope = 1,  linetype = "dashed", color = "gray", size = 0.5) +
    geom_abline(slope = (2/1)*(1/3),  linetype = "dashed", color = "indianred1", size = 0.5) +
    geom_abline(slope = (2/1)*(2/3),  linetype = "dashed", color = "cornflowerblue", size = 0.5) +
    
    geom_point(aes(color = Sex), alpha = 0.75, size = 0.3)+
  
    geom_smooth(method = "lm", alpha = 0.3, size = 0.5,color="black")+
    geom_smooth(method = "lm", alpha = 0.3, size = 0.5, aes(color = Sex))+
    guides(color=FALSE)+
    scale_color_manual(values = c("cornflowerblue", "indianred1"))+
    scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
    labs(x = "Autosomal expected \ngenetic contrib. in 2013", y = "Z expected genetic contrib. in 2013")+
    plottheme + 
    ylim(0,0.04) 


plot_grid(A_contribs_vs_descendants, Z_contribs_vs_descendants, A_vs_Z_contribs,  ncol = 3, labels = "AUTO", align = 'hv',axis='tblr')
#ggsave('fig2_EGC_AZ_6x2.pdf', width = 12, height = 4, units ='in')




