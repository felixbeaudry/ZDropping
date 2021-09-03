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


## Plot of offspring sex ratio vs. Z contribution for each mom

# read in pedigree
pedigree <- read.table("working_files/pedigree.txt", header = TRUE, stringsAsFactors = FALSE)

# squish pedigree to get sex ratio of moms' offspring, then combine with indiv contribs
# get rid of entries where mom is unknown
indiv_contribs_offspring_sex_ratio <- filter(pedigree, Mom != "0") %>%
  # group by mom and sex
  group_by(Mom, Sex) %>%
  # get number of offspring of each sex for each mom
  dplyr::summarize(num_offspring = n()) %>%
  # remove unsexed offspring as we can't do anything with these
  filter(Sex != 0) %>%
  # make sons and daughters into separate columns so that we can work with them
  pivot_wider(id_cols = Mom, names_from = Sex, names_prefix = "num_offspring_", values_from = num_offspring) %>%
  # calculate sex ratio: prop male offspring
  mutate(prop_males = num_offspring_1/(num_offspring_1+num_offspring_2)) %>%
  # pull out just the moms' IDs and sex ratios
  select(Indiv = Mom, prop_males) %>%
  # combine with indiv contribs data
  left_join(indiv_contribs_prop_nestlings_2013, .)

# plot
ggplot(indiv_contribs_offspring_sex_ratio, aes(x = prop_males, y = Z_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Offspring sex ratio (proportion male)", y = "Z expected genetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) +
  ylim(c(0,0.025)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.33, 0.5, 0.66, 0.75, 1))

# auto for comparison
ggplot(indiv_contribs_offspring_sex_ratio, aes(x = prop_males, y = auto_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Offspring sex ratio (proportion male)", y = "Autosomal expected genetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) +
  ylim(c(0,0.025)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.33, 0.5, 0.66, 0.75, 1))

# What about males
indiv_contribs_offspring_sex_ratio_dads <- filter(pedigree, Dad != "0") %>%
  # group by dad and sex
  group_by(Dad, Sex) %>%
  # get number of offspring of each sex for each dad
  dplyr::summarize(num_offspring = n()) %>%
  # remove unsexed offspring as we can't do anything with these
  filter(Sex != 0) %>%
  # make sons and daughters into separate columns so that we can work with them
  pivot_wider(id_cols = Dad, names_from = Sex, names_prefix = "num_offspring_", values_from = num_offspring) %>%
  # calculate sex ratio: prop male offspring
  mutate(prop_males = num_offspring_1/(num_offspring_1+num_offspring_2)) %>%
  # pull out just the moms' IDs and sex ratios
  select(Indiv = Dad, prop_males) %>%
  # combine with indiv contribs data
  left_join(indiv_contribs_prop_nestlings_2013, .)

# plot Z
ggplot(indiv_contribs_offspring_sex_ratio_dads, aes(x = prop_males, y = Z_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Offspring sex ratio (proportion male)", y = "Z expected genetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) +
  ylim(c(0,0.025)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.33, 0.5, 0.66, 0.75, 1))

# plot auto
ggplot(indiv_contribs_offspring_sex_ratio_dads, aes(x = prop_males, y = auto_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Offspring sex ratio (proportion male)", y = "Autosomal expected genetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) +
  ylim(c(0,0.025)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.33, 0.5, 0.66, 0.75, 1))


## Plot of male contrib vs. female contrib for all breeding pairs

# squish pedigree to get breeding pairs, then combine with indiv contribs
indiv_contribs_dad_vs_mom <- distinct(pedigree, Dad, Mom) %>%
  filter(Dad != "0", Mom != "0") %>%
  left_join(select(indiv_contribs_prop_nestlings_2013, Dad = Indiv, Dad_Z_mean = Z_mean, Dad_auto_mean = auto_mean)) %>%
  left_join(select(indiv_contribs_prop_nestlings_2013, Mom = Indiv, Mom_Z_mean = Z_mean, Mom_auto_mean = auto_mean))

# plot Z contribs
ggplot(indiv_contribs_dad_vs_mom, aes(x = Dad_Z_mean, y = Mom_Z_mean)) +
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  geom_point(alpha = 0.75, size = 0.3)+
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  labs(x = "Male's Z expected genetic contrib. in 2013", y = "Female's Z expected genetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))

# plot auto contribs for contrast
ggplot(indiv_contribs_dad_vs_mom, aes(x = Dad_auto_mean, y = Mom_auto_mean)) +
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  geom_point(alpha = 0.75, size = 0.3)+
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  labs(x = "Male's autosomal expected genetic contrib. in 2013", y = "Female's autosomal expected genetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))

# how many different partners do you have
dad_partners <- group_by(indiv_contribs_dad_vs_mom, Dad) %>%
  dplyr::summarize(dad_partners_n = n())
mom_partners <- group_by(indiv_contribs_dad_vs_mom, Mom) %>%
  dplyr::summarize(mom_partners_n = n())
indiv_contribs_dad_vs_mom_count_partners <- left_join(indiv_contribs_dad_vs_mom, dad_partners) %>%
  left_join(mom_partners) %>%
  mutate(partners_difference = dad_partners_n - mom_partners_n)

# distribution of partner numbers
ggplot(indiv_contribs_dad_vs_mom_count_partners, aes(x=dad_partners_n)) +
  geom_histogram()
ggplot(indiv_contribs_dad_vs_mom_count_partners, aes(x=mom_partners_n)) +
  geom_histogram()
ggplot(indiv_contribs_dad_vs_mom_count_partners, aes(x=partners_difference)) +
  geom_histogram()
ggplot(indiv_contribs_dad_vs_mom_count_partners, aes(x=dad_partners_n, y=mom_partners_n)) +
  geom_point(alpha = 0.01)

# pull out monogamous pairs
monogamous_pairs <- filter(indiv_contribs_dad_vs_mom_count_partners, dad_partners_n == 1, mom_partners_n == 1)
# plot
ggplot(monogamous_pairs, aes(x = Dad_Z_mean, y = Mom_Z_mean)) +
geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
geom_point(alpha = 0.75, size = 0.3)+
geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
labs(x = "Male's Z expected genetic contrib. in 2013", y = "Female's Z expected genetic contrib. in 2013")+
plottheme+
theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))
#auto
ggplot(monogamous_pairs, aes(x = Dad_auto_mean, y = Mom_auto_mean)) +
geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
geom_point(alpha = 0.75, size = 0.3)+
geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
labs(x = "Male's Z expected genetic contrib. in 2013", y = "Female's Z expected genetic contrib. in 2013")+
plottheme+
theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))

# combine with sex ratio of offspring
monogamous_pairs_offspring_sex_ratio <- left_join(monogamous_pairs, select(indiv_contribs_offspring_sex_ratio, Mom = Indiv, prop_males))
# plot Z
ggplot(monogamous_pairs_offspring_sex_ratio, aes(x = Dad_Z_mean, y = Mom_Z_mean)) +
#geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
geom_point(aes(color = prop_males), alpha = 0.75, size = 1)+
geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
labs(x = "Male's Z expected genetic contrib. in 2013", y = "Female's Z expected genetic contrib. in 2013")+
plottheme+
theme(plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))
# plot auto
ggplot(monogamous_pairs_offspring_sex_ratio, aes(x = Dad_auto_mean, y = Mom_auto_mean)) +
#geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
geom_point(aes(color = prop_males), alpha = 0.75, size = 0.3)+
geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
labs(x = "Male's Z expected genetic contrib. in 2013", y = "Female's Z expected genetic contrib. in 2013")+
plottheme+
theme(plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))
# who are those 6 pairs on the Z
filter(monogamous_pairs_offspring_sex_ratio, (Dad_Z_mean - Mom_Z_mean)>0.001)


# male vs female contribs, now colored by difference in number of partners

# plot Z
ggplot(indiv_contribs_dad_vs_mom_count_partners, aes(x = Dad_Z_mean, y = Mom_Z_mean, color = partners_difference)) +
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  geom_point(alpha = 0.75, size = 1)+
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  labs(x = "Male's Z expected genetic contrib. in 2013", y = "Female's Z expected genetic contrib. in 2013")+
  plottheme+
  theme(plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) +
  scale_color_gradient2(low = "darkred", mid = "mediumpurple", high = "darkblue")
# plot auto
ggplot(indiv_contribs_dad_vs_mom_count_partners, aes(x = Dad_auto_mean, y = Mom_auto_mean, color = partners_difference)) +
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5)+
  geom_point(alpha = 0.75, size = 1)+
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  labs(x = "Male's autosomal expected genetic contrib. in 2013", y = "Female's autosomal expected genetic contrib. in 2013")+
  plottheme+
  theme(plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))+
  scale_color_gradient2(low = "darkred", mid = "mediumpurple", high = "darkblue")

# male-female contribs vs difference in number of partners
indiv_contribs_dad_vs_mom_count_partners_exp_gen_contrib_diff <- mutate(indiv_contribs_dad_vs_mom_count_partners, Z_mean_diff = Dad_Z_mean - Mom_Z_mean, auto_mean_diff = Dad_auto_mean - Mom_auto_mean)
# Z
ggplot(indiv_contribs_dad_vs_mom_count_partners_exp_gen_contrib_diff, aes(x = partners_difference, y=Z_mean_diff)) +
  geom_point()+
  scale_y_continuous(limits = c(-0.02,0.03)) +
  geom_smooth(method = "lm")
# auto
ggplot(indiv_contribs_dad_vs_mom_count_partners_exp_gen_contrib_diff, aes(x = partners_difference, y=auto_mean_diff)) +
  geom_point() +
  scale_y_continuous(limits = c(-0.02,0.03)) +
  geom_smooth(method = "lm")

# fixing some colors
indiv_contribs_dad_vs_mom_count_partners$partners_difference_fct <- as.factor(indiv_contribs_dad_vs_mom_count_partners$partners_difference)

## Felix, have fun!