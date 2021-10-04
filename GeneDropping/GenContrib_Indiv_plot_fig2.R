#03 October 2021
#Rose Driscoll and Felix Beaudry
#script to plot the genetic & genealogical contributions of male and female breeders for autosomes and Z
#simplified version with just the plotting commands

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(kinship2)
library(cowplot)

#plot theme
plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3),
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"),
                    plot.background = element_rect(fill = "white"),
                    axis.text.x = element_text(size=6), axis.text.y = element_text(size=6),
                    axis.title = element_text(size=7), plot.title = element_text(size=8),
                    legend.position="right", legend.text = element_text(size=6),
                    legend.title = element_text(size=7), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(1,"cm"))


## Read in data

# Read in 1990-2013 breeders data
breeders_926 <- read.table("working_files/926_breeders.txt", header = TRUE, stringsAsFactors = FALSE)

# Read in pedigree
pedigree <- read.table("working_files/pedigree.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

# Read in all of the data tables including breeders data
load("all_tables_v2.rdata")

## Genealogical contributions for all breeders

# Make a list of the descendants of every bird
descendants <- NULL
for (i in breeders_926$Indiv) {
  generation = 2
  thisGen <- filter(pedigree, Mom==i | Dad==i) 
  repeat {
    if (nrow(thisGen)==0){
      break
    }
    thisGenRow <- select(thisGen, Indiv) %>%
      transmute(indiv_id = i, descendant_id = Indiv, relationship = generation)
    descendants <- rbind(descendants, thisGenRow)
    generation = generation + 1
    nextGen <- filter(pedigree, Mom %in% thisGen$Indiv | Dad %in% thisGen$Indiv)
    thisGen <- nextGen
  }
}

# get indiv data including natal years
indiv_data<-read.table('working_files/IndivDataUSFWS.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
nestling_data<-select(indiv_data, Indiv, NatalYear)

# add natal year info to descendants dataset
descendants_years <- left_join(descendants, nestling_data, by = c("descendant_id" = "Indiv"))

# how many nestlings are there each year?
num_nestlings_yearly <- group_by(nestling_data, NatalYear) %>%
  dplyr::summarize(num_nestlings = n())


# Get 2013 genetic contributions (autosomes and Z) for each individual

# Make a table to save individuals' 2013 contributions in
indiv_contribs_A_Z_2013<-NULL

# Only running this for breeders that actually have kids
breeders_926$Kids <- breeders_926$Indiv %in% c(pedigree$Dad, pedigree$Mom)

# Go through all of the individuals
for (i in breeders_926[breeders_926$Kids,'Indiv']) {
  
  #get autosomal data for individual i
  obsIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.A.1.drop.data.txt',sep=''),header=TRUE)
  simIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.A.1.drop.sim.txt',sep=''),header=TRUE)
  
  #make sure number of sims adds up to 1000000
  ddply(simIndiv[simIndiv$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))
  for (y in c(0:23)){
    if (sum(simIndiv[simIndiv$allele==2 & simIndiv$cohort==y,'all_alleles_count'])<1000000) {
      simIndiv<-rbind(simIndiv,cbind(cohort=y,allele=2,allele_count=0,all_alleles_count=1000000-sum(simIndiv[simIndiv$allele==2 & simIndiv$cohort==y,'all_alleles_count'])))
    }
  }
  ddply(simIndiv[simIndiv$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))
  
  #add year info
  simIndiv$Year<-simIndiv$cohort + 1990
  
  #subset data to only include 2013 and allele 2
  simIndiv<-simIndiv[simIndiv$allele==2 & simIndiv$Year == 2013,]
  
  #get total number of alleles sampled in 2013
  simIndiv$totAlleles<-laply(c(1:nrow(simIndiv)), function(x) unique(obsIndiv[obsIndiv$cohort_year==simIndiv[x,'Year'],'allele_count']))
  
  #summarize data
  simIndivSum<-ddply(simIndiv,.(Year),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), 
                     q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), 
                     q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))
  
  #save mean autosomal contribution for 2013
  indivauto <- mutate(simIndivSum, Indiv = i) %>%
    select(Indiv, auto_mean = mean)
  
  
  #get Z data for individual i
  obsIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.Z.1.drop.data.txt',sep=''),header=TRUE)
  simIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.Z.1.drop.sim.txt',sep=''),header=TRUE)
  
  #make sure number of sims adds up to 1000000
  ddply(simIndiv[simIndiv$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))
  for (y in c(0:23)){
    if (sum(simIndiv[simIndiv$allele==2 & simIndiv$cohort==y,'all_alleles_count'])<1000000) {
      simIndiv<-rbind(simIndiv,cbind(cohort=y,allele=2,allele_count=0,all_alleles_count=1000000-sum(simIndiv[simIndiv$allele==2 & simIndiv$cohort==y,'all_alleles_count'])))
    }
  }
  ddply(simIndiv[simIndiv$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))
  
  #add year info
  simIndiv$Year<-simIndiv$cohort + 1990
  
  #subset data to only include 2013 and allele 2
  simIndiv<-simIndiv[simIndiv$allele==2 & simIndiv$Year == 2013,]
  
  #get total number of alleles sampled in 2013
  simIndiv$totAlleles<-laply(c(1:nrow(simIndiv)), function(x) unique(obsIndiv[obsIndiv$cohort_year==simIndiv[x,'Year'],'allele_count']))
  
  #summarize data
  simIndivSum<-ddply(simIndiv,.(Year),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), 
                     q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), 
                     q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))
  
  #save mean Z contribution for 2013
  indivrow <- mutate(simIndivSum, Indiv = i) %>%
    select(Indiv, Z_mean = mean) %>%
    right_join(indivauto)
  
  # add individual to table
  indiv_contribs_A_Z_2013<-rbind(indiv_contribs_A_Z_2013, indivrow)
  
}

# Add back in breeders with no kids, with 0 contributions
indiv_contribs_A_Z_2013 <- left_join(select(breeders_926, Indiv), indiv_contribs_A_Z_2013) %>%
  mutate(Z_mean = ifelse(is.na(Z_mean), 0, Z_mean), auto_mean = ifelse(is.na(auto_mean), 0, auto_mean))

# add sex info from pedigree
indiv_contribs_A_Z_2013_sex <- left_join(indiv_contribs_A_Z_2013, select(pedigree, Indiv, Sex))
indiv_contribs_A_Z_2013_sex$Sex <- as.factor(indiv_contribs_A_Z_2013_sex$Sex)

# get total number of nestlings in 2013
num_2013_nestlings = filter(num_nestlings_yearly, NatalYear == 2013)$num_nestlings
# filter for descendants born in 2013 and summarize number of 2013 descendants for each founder
# then calculate what proportion of 2013 nestlings that is
descendants_2013_prop_nestlings <- filter(descendants_years, NatalYear == 2013) %>%
  group_by(indiv_id) %>%
  dplyr::summarize(num_descendants_2013 = n()) %>%
  left_join(dplyr::select(breeders_926, Indiv), ., by = c("Indiv" = "indiv_id")) %>%
  select(indiv_id = Indiv, num_descendants_2013) %>%
  mutate(num_descendants_2013 = ifelse(is.na(num_descendants_2013), 0, num_descendants_2013), total_2013_nestlings = num_2013_nestlings, prop_2013_nestlings = num_descendants_2013/total_2013_nestlings)

# join this data to the individual contributions data
indiv_contribs_prop_nestlings_2013 <- left_join(indiv_contribs_A_Z_2013_sex, descendants_2013_prop_nestlings, by = c("Indiv" = "indiv_id"))



####Pop Level distributions####
#genetic (A & Z) vs. genealogical contributions in 2013 for all 926 breeders

#count of individuals with zero descendants in 2013
sum(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==1] == 0)
sum(indiv_contribs_prop_nestlings_2013$prop_2013_nestlings[indiv_contribs_prop_nestlings_2013$Sex ==2] == 0)

## genealogical contributions
  #range and mean genealogical contributions for males with >0 descendants in 2013
  #range and mean genealogical contributions for females with >0 descendants in 2013
###expected genetic contributions
  #range and mean expected contributions for males with >0 descendants in 2013 - autosomal & Z
  #range and mean expected contributions for females with >0 descendants in 2013 - autosomal & Z
filter(indiv_contribs_prop_nestlings_2013, prop_2013_nestlings != 0) %>%
  group_by(Sex) %>%
  dplyr::summarize(geno_min = min(prop_2013_nestlings), geno_max = max(prop_2013_nestlings), geno_mean = mean(prop_2013_nestlings),
                   auto_min = min(auto_mean), auto_max = max(auto_mean), auto_mean = mean(auto_mean),
                   Z_min = min(Z_mean), Z_max = max(Z_mean), Z_mean = mean(Z_mean))


####descents vs genes####
#genealogical vs autosomal expected model
A_FM_mod_lm <- lm(auto_mean ~  prop_2013_nestlings*Sex  , data=indiv_contribs_prop_nestlings_2013)  
summary(A_FM_mod_lm)



#genealogical vs Z expected model
Z_FM_mod_lm <- lm(Z_mean ~  prop_2013_nestlings*Sex  , data=indiv_contribs_prop_nestlings_2013)  
summary(Z_FM_mod_lm)


#plot
(A_contribs_vs_descendants <- 
  ggplot(indiv_contribs_prop_nestlings_2013, aes(x = prop_2013_nestlings, y = auto_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5,formula=y~0+x)+
  geom_abline(slope = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Genealogical contribution in 2013", y = "Autosomal expected \ngenetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) + 
  ylim(0,0.04))


(Z_contribs_vs_descendants <- 
  ggplot(indiv_contribs_prop_nestlings_2013, aes(x = prop_2013_nestlings, y = Z_mean, color = Sex, fill = Sex)) +
  geom_point(alpha = 0.75, size = 0.3)+
  geom_smooth(method = "lm", alpha = 0.3, size = 0.5,formula=y~0+x)+
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("cornflowerblue", "indianred1"))+
  scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
  labs(x = "Genealogical contribution in 2013", y = "Z expected genetic contrib. in 2013")+
  plottheme+
  theme(legend.position = "none", plot.margin = unit(c(15,5.5,5.5,5.5), "pt")) + 
  ylim(0,0.04))

(A_vs_Z_contribs <- 
  ggplot(indiv_contribs_A_Z_2013_sex, aes(x = auto_mean, y = Z_mean)) +
    geom_abline(slope = 1,  linetype = "dashed", color = "gray", size = 0.5) +
    geom_abline(slope = (2/1)*(1/3),  linetype = "dashed", color = "indianred1", size = 0.5) +
    geom_abline(slope = (2/1)*(2/3),  linetype = "dashed", color = "cornflowerblue", size = 0.5) +
    geom_point(aes(color = Sex), alpha = 0.75, size = 0.3)+
    geom_smooth(method = "lm", alpha = 0.3, size = 0.5,color="black",formula=y~0+x)+
    geom_smooth(method = "lm", alpha = 0.3, size = 0.5, aes(color = Sex),formula=y~0+x)+
    guides(color=FALSE)+
    scale_color_manual(values = c("cornflowerblue", "indianred1"))+
    scale_fill_manual(values = c("cornflowerblue", "indianred1"))+
    labs(x = "Autosomal expected \ngenetic contrib. in 2013", y = "Z expected genetic contrib. in 2013")+
    plottheme + 
    ylim(0,0.04))


plot_grid(A_contribs_vs_descendants, Z_contribs_vs_descendants, A_vs_Z_contribs,  ncol = 3, labels = "AUTO", align = 'hv',axis='tblr')
#ggsave('fig2_EGC_AZ_6-5x2-26.pdf', width = 6.5, height = 2.26, units ='in')

#AZ

AZ_mod <- lm(Z_mean ~  auto_mean + 0  , data=indiv_contribs_A_Z_2013_sex)  
summary(AZ_mod)  

AZ_mods_sex <- lm(Z_mean ~  auto_mean*Sex + 0 , data=indiv_contribs_A_Z_2013_sex)  
summary(AZ_mods_sex)  

## Plot of male contrib vs. female contrib for breeding pairs amongst the 926 breeders

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
  # replace any NAs with 0
  mutate(num_offspring_1 = ifelse(is.na(num_offspring_1), 0, num_offspring_1), num_offspring_2 = ifelse(is.na(num_offspring_2), 0, num_offspring_2)) %>%
  # calculate sex ratio: prop male offspring
  mutate(prop_males = num_offspring_1/(num_offspring_1+num_offspring_2)) %>%
  # pull out just the moms' IDs and sex ratios
  select(Indiv = Mom, prop_males) %>%
  # combine with indiv contribs data; this will cut the list down to just moms in the 926 breeders list who have at least one offspring of known sex
  inner_join(indiv_contribs_prop_nestlings_2013)



# squish pedigree to get breeding pairs, then combine with indiv contribs
indiv_contribs_dad_vs_mom <- distinct(pedigree, Dad, Mom) %>%
  filter(Dad != "0", Mom != "0") %>%
  left_join(select(indiv_contribs_prop_nestlings_2013, Dad = Indiv, Dad_Z_mean = Z_mean, Dad_auto_mean = auto_mean)) %>%
  left_join(select(indiv_contribs_prop_nestlings_2013, Mom = Indiv, Mom_Z_mean = Z_mean, Mom_auto_mean = auto_mean))
# This table includes some pairs with no contribs data as they are not among the 926 breeders, so pull out just complete cases
indiv_contribs_dad_vs_mom <- indiv_contribs_dad_vs_mom[complete.cases(indiv_contribs_dad_vs_mom),]

# how many different partners do each of these birds have
dad_partners <- group_by(indiv_contribs_dad_vs_mom, Dad) %>%
  dplyr::summarize(dad_partners_n = n())
mom_partners <- group_by(indiv_contribs_dad_vs_mom, Mom) %>%
  dplyr::summarize(mom_partners_n = n())
indiv_contribs_dad_vs_mom_count_partners <- left_join(indiv_contribs_dad_vs_mom, dad_partners) %>%
  left_join(mom_partners) %>%
  mutate(partners_difference = dad_partners_n - mom_partners_n)

# pull out monogamous pairs
monogamous_pairs <- filter(indiv_contribs_dad_vs_mom_count_partners, dad_partners_n == 1, mom_partners_n == 1)

# combine with sex ratio of offspring
monogamous_pairs_offspring_sex_ratio <- inner_join(monogamous_pairs, select(indiv_contribs_offspring_sex_ratio, Mom = Indiv, prop_males))
# this reduces the monogamous pairs list to only those who have at least one offspring of known sex

# fixing some colors
indiv_contribs_dad_vs_mom_count_partners$partners_difference_fct <- as.factor(indiv_contribs_dad_vs_mom_count_partners$partners_difference)
legend_size=0.5
library(stringr)

pdf("fig_Sx_zeroContribs.pdf",width=5.5,height=4)

plot_grid(
  # offspring sex ratio; plot auto
  ggplot(monogamous_pairs_offspring_sex_ratio, aes(x = Dad_auto_mean, y = Mom_auto_mean)) +
    geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
    geom_point(aes(color = prop_males),  size = 1.2)+
    labs(x = "Male's autosomal expected genetic contrib. in 2013", y =str_wrap( "Female's autosomal expected genetic contrib. in 2013", width = 25),color="Progeny\nSex Ratio\n(Male Prop.)")+
    plottheme+
    theme(plot.margin = unit(c(15,5.5,5.5,5.5), "pt"),legend.key.size = unit(legend_size, 'cm'))+ scale_color_gradient(high = "cornflowerblue",  low = "indianred1"),
  
  # offspring sex ratio; plot Z
  ggplot(monogamous_pairs_offspring_sex_ratio, aes(x = Dad_Z_mean, y = Mom_Z_mean)) +
    geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
    geom_point(aes(color = prop_males),  size = 1.2)+
    labs(x = "Male's Z expected genetic contrib. in 2013", y = str_wrap("Female's Z expected genetic contrib. in 2013", width = 25),color="Progeny\nSex Ratio\n(Male Prop.)")+
    plottheme+
    theme(plot.margin = unit(c(15,5.5,5.5,5.5), "pt"),legend.key.size = unit(legend_size, 'cm')) + scale_color_gradient(high = "cornflowerblue",  low = "indianred1"),

  # partner diff; plot auto
  ggplot(indiv_contribs_dad_vs_mom_count_partners, aes(x = Dad_auto_mean, y = Mom_auto_mean, color = partners_difference)) +
    geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
    geom_point(size = 1.2)+
    labs(x = "Male's autosomal expected genetic contrib. in 2013", y = str_wrap("Female's autosomal expected genetic contrib. in 2013", width = 25),color="Partner Diff.\n(Male Bias)")+
    plottheme+
    theme(plot.margin = unit(c(15,5.5,5.5,5.5), "pt"),legend.key.size = unit(legend_size, 'cm'))+
    scale_color_gradient2(high = "cornflowerblue",  mid = "gray", low = "indianred1"),
  
  # partner diff; plot Z
  ggplot(indiv_contribs_dad_vs_mom_count_partners, aes(x = Dad_Z_mean, y = Mom_Z_mean, color = partners_difference)) +
    geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
    geom_point(size = 1.2)+
    labs(x = "Male's Z expected genetic contrib. in 2013", y = str_wrap("Female's Z expected genetic contrib. in 2013", width = 25),color="Partner Diff.\n(Male Bias)")+
    plottheme+
    theme(plot.margin = unit(c(15,5.5,5.5,5.5), "pt"),legend.key.size = unit(legend_size, 'cm')) +
    scale_color_gradient2(high = "cornflowerblue",  mid = "gray", low = "indianred1") ,
          
 
          ncol = 2, labels = "AUTO", align = 'hv',axis='tblr')

dev.off()
