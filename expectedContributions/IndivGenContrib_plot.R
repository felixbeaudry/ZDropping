#29 May 2019
#Rose Driscoll
#script to plot the genetic & genealogical contributions of male and female breeders for autosomes and Z
#simplified version with just the plotting commands

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(kinship2)
library(cowplot)

# load Rdata file 
setwd('~/Google Drive/Research/Data2/fsj/Zdropping')

load("plotIndivGenContrib_tidy_20190529.Rdata")
# code to produce all of these tables from raw data is in plotIndivGenContrib_tidy_20190425.R

# Nancy's plot theme

font = 12

plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3), 
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"), 
                    plot.background = element_rect(fill = "white"),  
                    axis.text.x = element_text(size=font), 
                    axis.text.y = element_text(size=font), 
                    axis.title = element_text(size=font+1), 
                    plot.title = element_text(size=font+2), 
                    legend.position="right", legend.text = element_text(size=7),
                    legend.title = element_text(size=font+2), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(1,"cm"))


####Pop Level distrbutions####
#genetic (A & Z) vs. genealogical contributions in 2013 for all 926 breeders

# Part A: genealogical vs. autosomal genetic contributions

# Plot the number of descendants against autosomal genetic contributions
ggplot(indiv_contribs_prop_nestlings_2013, aes(x = prop_2013_nestlings, color = Sex, fill = Sex)) +
  geom_histogram(position="dodge")+
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Genealogical contribution in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))

length(is.na(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2]))
length(is.na(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1]))

mean(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2]))#females
mean(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1]))#males

min(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2]))
min(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1]))

max(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2]))
max(na.omit(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1]))


ggplot(indiv_contribs_prop_nestlings_2013, aes(x = auto_mean, color = Sex, fill = Sex)) +
  geom_histogram(position="dodge")+
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Exp. A. gen. contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) 

mean(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==2]))#females
mean(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==1]))#males

min(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==2 & !is.na(indiv_contribs_prop_nestlings_2013$num_descendants_2013)]))
min(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==1 & !is.na(indiv_contribs_prop_nestlings_2013$num_descendants_2013)]))

max(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==2]))
max(na.omit(indiv_contribs_prop_nestlings_2013$auto_mean[indiv_contribs_prop_nestlings_2013$Sex ==1]))



ggplot(indiv_contribs_prop_nestlings_2013, aes(x = Z_mean, color = Sex, fill = Sex)) +
  geom_histogram(position="dodge")+
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Exp. Z. gen. contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) 



mean(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==2]))#females
mean(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==1]))#males

min(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==2 & !is.na(indiv_contribs_prop_nestlings_2013$num_descendants_2013)]))
min(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==1 & !is.na(indiv_contribs_prop_nestlings_2013$num_descendants_2013) ]))

max(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==2]))
max(na.omit(indiv_contribs_prop_nestlings_2013$Z_mean[indiv_contribs_prop_nestlings_2013$Sex ==1]))

####descents vs genes####
A_contribs_vs_descendants <- 
  ggplot(indiv_contribs_prop_nestlings_2013, aes(x = prop_2013_nestlings, y = auto_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.25, size = 0.3)+
  geom_abline(slope = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Genealogical contribution in 2013", y = "Exp. A. gen. contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) + ylim(0,0.04)

#A_FM_mod_lm <- lm(prop_2013_nestlings ~  0 + auto_mean  , data=indiv_contribs_A_Z_2013_sex[indiv_contribs_A_Z_2013_sex$Sex == 1,])  
A_FM_mod_lm <- lm(prop_2013_nestlings ~  auto_mean*Sex  , data=indiv_contribs_prop_nestlings_2013)  
summary(A_FM_mod_lm)

Z_contribs_vs_descendants <- 
  ggplot(indiv_contribs_prop_nestlings_2013, aes(x = prop_2013_nestlings, y = Z_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.25, size = 0.3)+
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Genealogical contribution in 2013", y = "Exp. Z gen. contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) + ylim(0,0.04)

Z_FM_mod_lm <- lm(prop_2013_nestlings ~  Z_mean*Sex  , data=indiv_contribs_prop_nestlings_2013)  
summary(Z_FM_mod_lm)

# Part C: autosomal genetic vs. Z genetic contributions

A_vs_Z_contribs <- 
  ggplot(indiv_contribs_A_Z_2013_sex, aes(x = auto_mean, y = Z_mean)) +
    geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
    geom_abline(slope = (2/1)*(1/3),  linetype = "dashed", color = "indianred1") +
    geom_abline(slope = (2/1)*(2/3),  linetype = "dashed", color = "cornflowerblue") +
    
    geom_point(aes(color = Sex), alpha = 0.75, size = 0.3)+
  
    geom_smooth(method = "lm", alpha = 0.25, size = 0.3,color="black")+
    geom_smooth(method = "lm", alpha = 0.25, size = 0.3, aes(color = Sex))+
    guides(color=FALSE)+
    scale_color_manual(values = c("cornflowerblue", "indianred1"))+
    scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
    labs(x = "Exp. A. genetic contrib. in 2013", y = "Exp. Z genetic contrib. in 2013")+
    plottheme + 
    ylim(0,0.04) 



#models
  


AZ_mod <- lm(Z_mean ~  auto_mean  , data=indiv_contribs_A_Z_2013_sex)  
summary(AZ_mod_m)  

AZ_mods_sex <- lm(Z_mean ~  auto_mean*Sex  , data=indiv_contribs_A_Z_2013_sex)  
summary(AZ_mods_sex)  


AZ_mod_m <- lm(Z_mean ~  auto_mean  , data=indiv_contribs_A_Z_2013_sex[indiv_contribs_A_Z_2013_sex$Sex == 1,])  
AZ_mod_f <- lm(Z_mean ~  auto_mean  , data=indiv_contribs_A_Z_2013_sex[indiv_contribs_A_Z_2013_sex$Sex == 2,])  
summary(AZ_mod_m)  
summary(AZ_mod_f)  


AZ_mod <- lm(auto_mean ~ Z_mean , data=indiv_contribs_A_Z_2013_sex)  
summary(AZ_mod)  

AZ_mod <- lm(Z_mean ~  auto_mean, data=indiv_contribs_A_Z_2013_sex)  
summary(AZ_mod)  




plot_grid(A_contribs_vs_descendants, Z_contribs_vs_descendants, A_vs_Z_contribs,  ncol = 3, labels = "AUTO", align = 'hv',axis='tblr')
ggsave('fig2_EGC_AZ_6x2.pdf', width = 12, height = 4, units ='in')




