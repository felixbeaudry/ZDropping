#script to plot the genetic contribution of male and female immigrant cohorts to the 2013 birth cohort
#vs the number of individuals in each cohort
#for immigrant cohorts that were all dead by 2013
#Rose Driscoll + some code from Nancy Chen
#7 April 2021

## Setup
setwd('~/Google Drive/Research/Data2/fsj/')

library(plyr)
library(ggplot2)
library(dplyr)

# plot theme
plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3), 
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"), 
                    plot.background = element_rect(fill = "white"),  
                    axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), 
                    axis.title = element_text(size=7), plot.title = element_text(size=8), 
                    legend.position="right", legend.text = element_text(size=7),
                    legend.title = element_text(size=8), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(1,"cm"))


## Immigrant cohorts

# First get a list of imms
# IndivDataUSFWS.txt says who is from which immigrant cohort
indivdata <- read.table("ZDropping_all/IndivDataUSFWS.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Just need the imms
immdata <- filter(indivdata, !is.na(ImmCohort)) %>%
  select(Indiv,ImmCohort)

# Now we need specifically imms who are breeders
# Find this out from the pedigree
ped <- read.table("ZDropping_all/FSJpedgeno_Zsexlinked.ped", header = FALSE, sep = " ", stringsAsFactors = FALSE)
pedigree <- ped[,1:6]
colnames(pedigree) <- c("Fam", "Indiv", "Dad", "Mom", "Sex", "Pheno")
# Which imms are breeders (find out by checking if they're someone's mom or dad)
immbreeders <- filter(immdata, Indiv %in% c(pedigree$Dad, pedigree$Mom))

# all.obs has observations of individuals
load("databases/all_tables_v2.rdata")
# find the last observation of each indiv
immbreeders <- filter(all.obs, USFWS %in% immbreeders$Indiv) %>%
  group_by(USFWS) %>%
  summarize(last_obs = max(Year)) %>%
  left_join(immbreeders, ., by = c("Indiv" = "USFWS"))

# How long does it take for an imm cohort to all be dead?
group_by(immbreeders, ImmCohort) %>%
  summarize(last_obs_of_cohort = max(last_obs)) %>%
  mutate(time_to_last_obs = last_obs_of_cohort-ImmCohort) %>%
  filter(ImmCohort>1999)
#ImmCohort last_obs_of_cohort time_to_last_obs
#1991                    2004               13
#1992                    2007               15
#1993                    2008               15
#1994                    2007               13
#1995                    2005               10
#1996                    2010               14
#1997                    2010               13
#1998                    2014               16
#1999                    2009               10
#2000                    2014               14
#2001                    2008                7
#2002                    2014               12
#2003                    2011                8
#2004                    2015               11
#2005                    2013                8
#2006                    2014                8
#2007                    2014                7
#2008                    2014                6
#2009                    2014                5
#2010                    2014                4
#2011                    2015                4
#2012                    2014                2
# Can't have anything where the last obs is 2014...
# Need to allow 15 years for them all to be dead due to 1992 & 1993
# So we can only go to 1998, but 1998 is omitted because the last obs is in 2014
# So we can only go to 1997

# Cohorts to use
imm_cohorts_to_use <- c(1991:1997)

# How many immigrants are in each cohort?
num_MF_imms_per_cohort <- left_join(immbreeders, pedigree) %>%
  group_by(ImmCohort, Sex) %>%
  summarize(num_imms = n()) %>% 
  filter(ImmCohort < 1998)



## Expected genetic contributions of immigrant cohorts to the birth cohort 15 years later

# Code for this section borrowed from plotImmGenContrib_malevsfemale_ZvsA.R

# Auto

#get data
obsImmYrA<-read.table('ZDropping_all/ImmContribYearly_malevsfemale_A.1.drop.data.txt',header=TRUE)
simImmYrA<-read.table('ZDropping_all/ImmContribYearly_malevsfemale_A.1.drop.sim.txt',header=TRUE)

#add year
simImmYrA$Year<-simImmYrA$cohort + 1990

#males

#subset alleles that count (males are evens) and truncate at 2013
simImmYrAM<-simImmYrA[simImmYrA$Year<2014 & simImmYrA$allele %in% seq(2,44,by=2),]

#make sure num sims adds up to 1000000
ddply(simImmYrAM,.(Year,allele),summarize, sum=sum(all_alleles_count))

for (yr in c(1990:2013)){
  for(a in unique(simImmYrAM[simImmYrAM$Year==yr,'allele'])) {
    if (sum(simImmYrAM[simImmYrAM$Year==yr & simImmYrAM$allele==a,'all_alleles_count'])<1000000) {
      simImmYrAM<-rbind(simImmYrAM,cbind(cohort=yr-1990,allele=a,allele_count=0,all_alleles_count=1000000-sum(simImmYrAM[simImmYrAM$Year==yr & simImmYrAM$allele==a,'all_alleles_count']),Year=yr))
    }		
  }
}

#get total number of alleles sampled each year
simImmYrAM$totAlleles<-laply(c(1:length(simImmYrAM$cohort)), function(x) unique(obsImmYrA[obsImmYrA$cohort_year==simImmYrAM[x,'Year'],'allele_count']))

#get mean 
simImmYrAAvgM<-ddply(simImmYrAM,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count,all_alleles_count))))

simsumImmYrAM<-ddply(simImmYrAM,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))
simsumImmYrAM<-simsumImmYrAM[simsumImmYrAM$allele != 1,]

#females

#subset alleles that count (females are odds) and truncate at 2013
simImmYrAF<-simImmYrA[simImmYrA$Year<2014 & simImmYrA$allele %in% seq(3,45,by=2),]

#make sure num sims adds up to 1000000
ddply(simImmYrAF,.(Year,allele),summarize, sum=sum(all_alleles_count))

for (yr in c(1990:2013)){
  for(a in unique(simImmYrAF[simImmYrAF$Year==yr,'allele'])) {
    if (sum(simImmYrAF[simImmYrAF$Year==yr & simImmYrAF$allele==a,'all_alleles_count'])<1000000) {
      simImmYrAF<-rbind(simImmYrAF,cbind(cohort=yr-1990,allele=a,allele_count=0,all_alleles_count=1000000-sum(simImmYrAF[simImmYrAF$Year==yr & simImmYrAF$allele==a,'all_alleles_count']),Year=yr))
    }		
  }
}

#get total number of alleles sampled each year
simImmYrAF$totAlleles<-laply(c(1:length(simImmYrAF$cohort)), function(x) unique(obsImmYrA[obsImmYrA$cohort_year==simImmYrAF[x,'Year'],'allele_count']))

#get mean 
simImmYrAAvgF<-ddply(simImmYrAF,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count,all_alleles_count))))

simsumImmYrAF<-ddply(simImmYrAF,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))
simsumImmYrAF<-simsumImmYrAF[simsumImmYrAF$allele != 1,]

# Z

#get data
obsImmYrZ<-read.table('ZDropping_all/ImmContribYearly_malevsfemale_Z.1.drop.data.txt',header=TRUE)
simImmYrZ<-read.table('ZDropping_all/ImmContribYearly_malevsfemale_Z.1.drop.sim.txt',header=TRUE)

#add year
simImmYrZ$Year<-simImmYrZ$cohort + 1990

#males

#subset alleles that count (males are evens) and truncate at 2013
simImmYrZM<-simImmYrZ[simImmYrZ$Year<2014 & simImmYrZ$allele %in% seq(2,44,by=2),]

#make sure num sims adds up to 1000000
ddply(simImmYrZM,.(Year,allele),summarize, sum=sum(all_alleles_count))

for (yr in c(1990:2013)){
  for(a in unique(simImmYrZM[simImmYrZM$Year==yr,'allele'])) {
    if (sum(simImmYrZM[simImmYrZM$Year==yr & simImmYrZM$allele==a,'all_alleles_count'])<1000000) {
      simImmYrZM<-rbind(simImmYrZM,cbind(cohort=yr-1990,allele=a,allele_count=0,all_alleles_count=1000000-sum(simImmYrZM[simImmYrZM$Year==yr & simImmYrZM$allele==a,'all_alleles_count']),Year=yr))
    }		
  }
}

#get total number of alleles sampled each year
simImmYrZM$totAlleles<-laply(c(1:length(simImmYrZM$cohort)), function(x) unique(obsImmYrZ[obsImmYrZ$cohort_year==simImmYrZM[x,'Year'],'allele_count']))

#get mean 
simImmYrZAvgM<-ddply(simImmYrZM,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count,all_alleles_count))))

simsumImmYrZM<-ddply(simImmYrZM,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))
simsumImmYrZM<-simsumImmYrZM[simsumImmYrZM$allele != 1,]

#females

#subset alleles that count (females are odds) and truncate at 2013
simImmYrZF<-simImmYrZ[simImmYrZ$Year<2014 & simImmYrZ$allele %in% seq(3,45,by=2),]

#make sure num sims adds up to 1000000
ddply(simImmYrZF,.(Year,allele),summarize, sum=sum(all_alleles_count))

for (yr in c(1990:2013)){
  for(a in unique(simImmYrZF[simImmYrZF$Year==yr,'allele'])) {
    if (sum(simImmYrZF[simImmYrZF$Year==yr & simImmYrZF$allele==a,'all_alleles_count'])<1000000) {
      simImmYrZF<-rbind(simImmYrZF,cbind(cohort=yr-1990,allele=a,allele_count=0,all_alleles_count=1000000-sum(simImmYrZF[simImmYrZF$Year==yr & simImmYrZF$allele==a,'all_alleles_count']),Year=yr))
    }		
  }
}

#get total number of alleles sampled each year
simImmYrZF$totAlleles<-laply(c(1:length(simImmYrZF$cohort)), function(x) unique(obsImmYrZ[obsImmYrZ$cohort_year==simImmYrZF[x,'Year'],'allele_count']))

#get mean 
simImmYrZAvgF<-ddply(simImmYrZF,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count,all_alleles_count))))

simsumImmYrZF<-ddply(simImmYrZF,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))
simsumImmYrZF<-simsumImmYrZF[simsumImmYrZF$allele != 1,]


## Grab the data I need for the plot and add imm cohort info

# convert allele to cohort year
# males have allele 2 = 1991, allele 4 = 1992, allele 6 = 1993, etc
# females have allele 3 = 1991, allele 5 = 1992, allele 7 = 1993, etc
# so for males it is allele/2 + 1990
# and for females it is (allele-1)/2 + 1990
simsumImmYrAFplus15 <- mutate(simsumImmYrAF, ImmCohort = ((allele-1)/2) + 1990, Sex = 2) %>%
  filter(Year==ImmCohort+15, ImmCohort < 1998)
simsumImmYrAMplus15 <- mutate(simsumImmYrAM, ImmCohort = (allele/2) + 1990, Sex = 1) %>%
  filter(Year==ImmCohort+15, ImmCohort < 1998)
simsumImmYrZFplus15 <- mutate(simsumImmYrZF, ImmCohort = ((allele-1)/2) + 1990, Sex = 2) %>%
  filter(Year==ImmCohort+15, ImmCohort < 1998)
simsumImmYrZMplus15 <- mutate(simsumImmYrZM, ImmCohort = (allele/2) + 1990, Sex = 1) %>%
  filter(Year==ImmCohort+15, ImmCohort < 1998)

simsumImmYrAplus15 <- bind_rows(simsumImmYrAFplus15, simsumImmYrAMplus15) %>%
  left_join(num_MF_imms_per_cohort, by = c("ImmCohort", "Sex")) %>%
  mutate(Sex = ifelse(Sex == 1, "M", "F"))
simsumImmYrZplus15 <- bind_rows(simsumImmYrZFplus15, simsumImmYrZMplus15) %>%
  left_join(num_MF_imms_per_cohort, by = c("ImmCohort", "Sex")) %>%
  mutate(Sex = ifelse(Sex == 1, "M", "F"))


## Plot!

# Num imms

# Auto
ggplot(simsumImmYrAplus15, aes(x=num_imms, y=mean, color=ImmCohort)) +
  geom_point(aes(shape = Sex)) +
  geom_linerange(aes(ymin=q1, ymax=q3)) +
  labs(x="Size of immigrant cohort, sex-specific \n(number of birds of the specified sex)", y="Expected autosomal genetic contribution to \nbirth cohort 15 years after immigration", title = "Autosomes")
ggplot(simsumImmYrAplus15, aes(x=num_imms, y=mean, color=Sex)) +
  geom_point(show.legend = FALSE) +
  geom_linerange(aes(ymin=q1, ymax=q3), show.legend = FALSE) +
  geom_smooth(method = "lm", aes(fill = Sex), alpha = 0.25) +
  scale_color_manual(values = c("indianred1", "cornflowerblue")) +
  scale_fill_manual(values = c("indianred1", "cornflowerblue")) +
  labs(x="Size of immigrant cohort, sex-specific \n(number of birds of the specified sex)", y="Expected autosomal genetic contribution to \nbirth cohort 15 years after immigration", title = "Autosomes") +
  plottheme +
  guides(color=guide_legend(override.aes=list(color=NA)))
  #https://stackoverflow.com/questions/35108443/ggplot2-make-legend-key-fill-transparent?rq=1
auto_lm <- lm(data = simsumImmYrAplus15, mean~num_imms+Sex)
summary(auto_lm)
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   -0.0111733  0.0173173  -0.645   0.5333  
#num_imms       0.0023911  0.0010354   2.309   0.0436 *
#SexM           0.0035813  0.0270021   0.133   0.8971  
#num_imms:SexM  0.0002614  0.0018748   0.139   0.8919  

# Z
ggplot(simsumImmYrZplus15, aes(x=num_imms, y=mean, color=ImmCohort)) +
  geom_point(aes(shape = Sex)) +
  geom_linerange(aes(ymin=q1, ymax=q3)) +
  labs(x="Size of immigrant cohort, sex-specific \n(number of birds of the specified sex)", y="Expected Z genetic contribution to \nbirth cohort 15 years after immigration", title = "Z")
ggplot(simsumImmYrZplus15, aes(x=num_imms, y=mean, color=Sex)) +
  geom_point(show.legend = FALSE) +
  geom_linerange(aes(ymin=q1, ymax=q3), show.legend = FALSE) +
  geom_smooth(method = "lm", aes(fill = Sex), alpha = 0.25) +
  scale_color_manual(values = c("indianred1", "cornflowerblue")) +
  scale_fill_manual(values = c("indianred1", "cornflowerblue")) +
  labs(x="Size of immigrant cohort, sex-specific \n(number of birds of the specified sex)", y="Expected Z genetic contribution to \nbirth cohort 15 years after immigration", title = "Z") +
  plottheme +
  guides(color=guide_legend(override.aes=list(color=NA)))
Z_lm <- lm(data = simsumImmYrZplus15, mean~Sex+num_imms)
summary(Z_lm)
0.019519/0.001616
#                Estimate Std. Error t value Pr(>|t|)   
#(Intercept)    0.001428   0.019102   0.075    0.942
#num_imms       0.001024   0.001142   0.897    0.391
#SexM          -0.007284   0.029785  -0.245    0.812
#num_imms:SexM  0.001940   0.002068   0.938    0.370
Z_M_lm <- lm(data = filter(simsumImmYrZplus15, Sex=="M"), mean~num_imms)
summary(Z_M_lm)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.005856   0.029711  -0.197    0.852
#num_imms     0.002964   0.002241   1.322    0.243
Z_F_lm <- lm(data = filter(simsumImmYrZplus15, Sex=="F"), mean~num_imms)
summary(Z_F_lm)
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 0.0014275  0.0106319   0.134    0.898
#num_imms    0.0010242  0.0006357   1.611    0.168
0.002964/0.0010242
#2.893966 # male slope is ~3x that of female slope, ??
ggplot(simsumImmYrZplus15, aes(x=num_imms, y=mean, color=Sex)) +
  geom_point(show.legend = FALSE) +
  geom_linerange(aes(ymin=q1, ymax=q3), show.legend = FALSE) +
  geom_smooth(method = "lm", aes(fill = Sex), alpha = 0.25) +
  scale_color_manual(values = c("indianred1", "cornflowerblue")) +
  scale_fill_manual(values = c("indianred1", "cornflowerblue")) +
  geom_abline(slope = 0.001616*(0.019519*100), intercept = -0.008040, linetype = "dashed") +
  geom_abline(slope = 0.001616, intercept = -0.008040, linetype = "longdash") +
  annotate("text", x=22, y=0.075, label = "Male slope = \n0.002964") +
  annotate("text", x=23, y=0.04, label = "1/2 of male slope") +
  labs(x="Size of immigrant cohort, sex-specific \n(number of birds of the specified sex)", y="Expected Z genetic contribution to \nbirth cohort 15 years after immigration", title = "Z") +
  plottheme +
  guides(color=guide_legend(override.aes=list(color=NA)))
ggplot(simsumImmYrZplus15, aes(x=num_imms, y=mean, color=Sex)) +
  geom_point(show.legend = FALSE) +
  geom_linerange(aes(ymin=q1, ymax=q3), show.legend = FALSE) +
  geom_smooth(method = "lm", aes(fill = Sex), alpha = 0.25) +
  scale_color_manual(values = c("indianred1", "cornflowerblue")) +
  scale_fill_manual(values = c("indianred1", "cornflowerblue")) +
  geom_abline(slope = 0.0010242, intercept = 0.0014275, linetype = "dashed") +
  geom_abline(slope = 0.0010242*2, intercept = 0.0014275, linetype = "dashed") +
  annotate("text", x=22, y=0.06, label = "Double \nfemale slope") +
  annotate("text", x=21, y=0.015, label = "Female slope = 0.0010242") +
  labs(x="Size of immigrant cohort, sex-specific \n(number of birds of the specified sex)", y="Expected Z genetic contribution to \nbirth cohort 15 years after immigration", title = "Z") +
  plottheme +
  guides(color=guide_legend(override.aes=list(color=NA)))


# save the data frames

save(simsumImmYrAplus15, simsumImmYrZplus15, file = "plotImmGenContrib_num_imms_vs_birth_cohort_contrib_20210407.rdata")
