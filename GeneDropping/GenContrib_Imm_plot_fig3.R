#script to plot the genetic contributions of male and female immigrants for autosomes and Z - fig 3
#Rose Driscoll and Felix Beaudry

library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(forcats)

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


####Part A: Number of imms (sexed & unsexed, breeding & nonbreeding) each year####

# Read in pedigree
pedigree <- read.table("working_files/pedigree.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
colnames(pedigree) <- c("Fam", "Indiv", "Dad", "Mom", "Sex", "Pheno")
pedigree$Indiv <- as.character(pedigree$Indiv)

# Read in immigrant cohort data
indivdata <- read.table("working_files/IndivData.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
indivdata$Indiv <- as.character(indivdata$Indiv)
# just need immigrant cohorts
immdata <- dplyr::select(indivdata, Indiv, ImmCohort)
# And join the immigrant cohorts to the pedigree
pedigree_immdata <- left_join(pedigree, immdata, by = "Indiv")




# Figure out which imms end up breeding
pedigree_immdata_breeder_nonbreeder_imms <- filter(pedigree_immdata, !is.na(ImmCohort)) %>%
  mutate(does_breed = ifelse(Indiv %in% pedigree_immdata$Dad, TRUE, 
                             ifelse(Indiv %in% pedigree_immdata$Mom, TRUE, FALSE)))
# Now count up how many imms that end up breeding or not come in each year
pedigree_immdata_breeder_nonbreeder_imms_count <- group_by(pedigree_immdata_breeder_nonbreeder_imms, ImmCohort, Sex, does_breed) %>%
  dplyr::reframe(num_imms = n()) %>%
  tidyr::complete(ImmCohort, Sex, does_breed, fill = list(num_imms = 0)) %>%
  filter(!is.na(ImmCohort))

pedigree_immdata_breeder_nonbreeder_imms_perc <- pedigree_immdata_breeder_nonbreeder_imms_count %>%
  group_by( ImmCohort,does_breed) %>% transmute(Sex, percent = num_imms/sum(num_imms))
# make males negative
pedigree_immdata_breeder_nonbreeder_imms_count_males_negative <- mutate(pedigree_immdata_breeder_nonbreeder_imms_count, num_imms_males_negative = ifelse(Sex==1, -num_imms, num_imms))

# combo plot
# plot sexed birds
(sexed_breed_nonbreed <- ggplot(filter(pedigree_immdata_breeder_nonbreeder_imms_count_males_negative, Sex != 0)) +
    geom_bar(aes(x = ImmCohort, y= num_imms_males_negative, fill = interaction(does_breed, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("lightskyblue1", "cornflowerblue", "mistyrose", "indianred1")) +
    labs(x = "", y = "Number of immigrants") +
    plottheme +
    scale_y_continuous(breaks = c(-30, -20, -10, 0, 10, 20, 30, 40, 50), labels = c(30, 20, 10, 0, 10, 20, 30, 40, 50)) +
    theme(legend.position = "none", plot.margin=unit(c(0.2,0.1,0,0.15),'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()))
# plot unsexed birds
(unsexed_breed_nonbreed <- ggplot(filter(pedigree_immdata_breeder_nonbreeder_imms_count_males_negative, Sex == 0)) +
    geom_bar(aes(x = ImmCohort, y= num_imms_males_negative, fill = interaction(does_breed, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("gray85")) +
    labs(x = "Year", y = "") +
    plottheme +
    scale_y_continuous(breaks = c(0, 10, 20)) +
    theme(legend.position = "none", plot.margin=unit(c(0.2,0.1,0,0.15),'cm')))
# plot to make a legend
for_legend <- mutate(pedigree_immdata_breeder_nonbreeder_imms_count_males_negative, fill_group = paste(does_breed, Sex, sep = "."))
for_legend$fill_group <- factor(for_legend$fill_group, levels = c("FALSE.2", "FALSE.1", "FALSE.0", "TRUE.2", "TRUE.1"))
(plot_for_legend <- ggplot(for_legend) +
    geom_bar(aes(x = ImmCohort, y= num_imms_males_negative, fill = fill_group), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("mistyrose", "lightskyblue1", "gray85", "indianred1", "cornflowerblue"), labels = c("Female\nnonbreeder", "Male\nnonbreeder", "Unsexed\nnonbreeder", "Female\nbreeder", "Male\nbreeder")) +
    labs(fill = "") +
    plottheme +
    theme(legend.key.size = unit(0.7,"cm"), legend.position = "top"))
sexed_unsexed_breed_nonbreed_legend <- get_legend(plot_for_legend)
# combine
sexed_unsexed_breed_nonbreed_plot <- plot_grid(sexed_breed_nonbreed, unsexed_breed_nonbreed, ncol = 1, nrow = 2, rel_heights = c(0.8, 0.25), align = "v")
sexed_unsexed_breed_nonbreed_plot_w_legend <- plot_grid(sexed_unsexed_breed_nonbreed_plot, sexed_unsexed_breed_nonbreed_legend, ncol = 2, nrow = 1, rel_widths = c(0.8, 0.25))

#linear model
pedigree_immdata_breeder_nonbreeder_imms_count$ImmCohort_zerod <- pedigree_immdata_breeder_nonbreeder_imms_count$ImmCohort - min(pedigree_immdata_breeder_nonbreeder_imms_count$ImmCohort)
immigrating_lm <- lm(data=pedigree_immdata_breeder_nonbreeder_imms_count, num_imms ~ Sex + ImmCohort_zerod)
summary(immigrating_lm)


#####Parts B & C: Autosomal & Z exp gen contribs of male & female imms#####

all_years_alleles <- data.frame(year = rep(c(1990:2013), each = 2), allele = rep(c(2,3), times = length(c(1990:2013))))

####Part B: Autosomal exp gen contribs of male & female imms####

# Read in data
obsImmA <- read.table('working_files/intermediate_files/ImmContribAll_prop_male_1_A.1.drop.data.txt', header=TRUE)
simImmA <- read.table('working_files/intermediate_files/ImmContribAll_prop_male_1_A.1.drop.sim.txt', header=TRUE)

# Add a column for year to the sim data
simImmA <- mutate(simImmA, year = cohort + 1990)

# Filter for alleles that count and truncate at 2013
simImmA <- filter(simImmA, allele > 1, year < 2014)

# Get total number of alleles each year from the obs data
totAlleles_dfA <- select(obsImmA, cohort, totAlleles = allele_count) %>%
  distinct()
simImmA <- left_join(simImmA, totAlleles_dfA, by = "cohort")

# Calculate means and quantiles
simsumImmA <- group_by(simImmA, year, allele) %>%
  dplyr::summarize(mean = mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), 
            q1 = quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), 
            q3 = quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))

# Make the expected genetic contributions plots start at 1990 by adding 0s for years where there is no data
simsumImmA_from1990 <- left_join(all_years_alleles, simsumImmA) %>%
  mutate(mean = ifelse(is.na(mean), 0, mean), q1 = ifelse(is.na(q1), 0, q1), q3 = ifelse(is.na(q3), 0, q3))

# plot
(A_imm_gen_contribs <- 
  ggplot(simsumImmA_from1990) + 
  geom_ribbon(aes(x = year, ymin = q1, ymax = q3, fill = as.factor(allele)), alpha=0.3) +
  geom_line(aes(x = year, y = mean, col = as.factor(allele)), size=0.5) + 
  labs(x = "Year", y = "Autosomal expected genetic contrib.") +
  scale_y_continuous(limits = c(0, 0.55), breaks = seq(0,0.5,by=0.1)) +
  scale_color_manual(values = c("cornflowerblue", "indianred1")) +
  scale_fill_manual(values = c("cornflowerblue", "indianred1")) +
  plottheme + 
  theme(legend.position='none',plot.margin=unit(c(0.2,0.1,0,0.15),'cm')))

#cumulative immigration proportion
simsumImmA_from1990$mean[simsumImmA_from1990$year == "2013" & simsumImmA_from1990$allele == 2] + 
  simsumImmA_from1990$mean[simsumImmA_from1990$year == "2013" & simsumImmA_from1990$allele == 3]  
  
  simsumImmA_from1990$mean[simsumImmA_from1990$year == "2013" & simsumImmA_from1990$allele == 3]/  
  (  simsumImmA_from1990$mean[simsumImmA_from1990$year == "2013" & simsumImmA_from1990$allele == 2] + 
    simsumImmA_from1990$mean[simsumImmA_from1990$year == "2013" & simsumImmA_from1990$allele == 3]  )
  

####Part C: Z exp gen contribs of male & female imms####

# Read in data
obsImmZ <- read.table('working_files/intermediate_files/ImmContribAll_prop_male_1_Z.1.drop.data.txt', header=TRUE)
simImmZ <- read.table('working_files/intermediate_files/ImmContribAll_prop_male_1_Z.1.drop.sim.txt', header=TRUE)

# Add a column for year to the sim data
simImmZ <- mutate(simImmZ, year = cohort + 1990)

# Filter for alleles that count and truncate at 2013
simImmZ <- filter(simImmZ, allele > 1, year < 2014)

# Get total number of alleles each year from the obs data
totAlleles_dfZ <- select(obsImmZ, cohort, totAlleles = allele_count) %>%
  distinct()
simImmZ <- left_join(simImmZ, totAlleles_dfZ, by = "cohort")

# Calculate means and quantiles
simsumImmZ <- group_by(simImmZ, year, allele) %>%
  dplyr::summarize(mean = mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), 
            q1 = quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), 
            q3 = quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))

# Make the expected genetic contributions plots start at 1990 by adding 0s for years where there is no data
simsumImmZ_from1990 <- left_join(all_years_alleles, simsumImmZ) %>%
  mutate(mean = ifelse(is.na(mean), 0, mean), q1 = ifelse(is.na(q1), 0, q1), q3 = ifelse(is.na(q3), 0, q3))

# plot
(Z_imm_gen_contribs <- 
  ggplot(simsumImmZ_from1990) + 
  geom_ribbon(aes(x = year, ymin = q1, ymax = q3, fill = as.factor(allele)), alpha=0.3) +
  geom_line(aes(x = year, y = mean, col = as.factor(allele)), size=0.5) + 
  labs(x = "Year", y = "Z expected genetic contrib.") +
  scale_y_continuous(limits = c(0, 0.55), breaks = seq(0,0.5,by=0.1)) +
  scale_color_manual(values = c("cornflowerblue", "indianred1")) +
  scale_fill_manual(values = c("cornflowerblue", "indianred1")) +
  plottheme + 
  theme(legend.position='none',plot.margin=unit(c(0.2,0.1,0,0.15),'cm')))


#cumulative immigration proportion
simsumImmZ_from1990$mean[simsumImmZ_from1990$year == "2013" & simsumImmZ_from1990$allele == 2] + 
  simsumImmZ_from1990$mean[simsumImmZ_from1990$year == "2013" & simsumImmZ_from1990$allele == 3]  

simsumImmZ_from1990$mean[simsumImmZ_from1990$year == "2013" & simsumImmZ_from1990$allele == 3]/  
  (  simsumImmZ_from1990$mean[simsumImmZ_from1990$year == "2013" & simsumImmZ_from1990$allele == 2] + 
       simsumImmZ_from1990$mean[simsumImmZ_from1990$year == "2013" & simsumImmZ_from1990$allele == 3]  )


#####Parts D & E: Autosomal & Z exp gen contribs 15 years after immigration vs imm cohort size#####

## Immigrant cohorts

# First get a list of imms
immdata_imms <- filter(indivdata, !is.na(ImmCohort)) %>%
  select(Indiv,ImmCohort)

# Now we need specifically imms who are breeders
# Find this out from the pedigree by checking if they're someone's mom or dad
immbreeders <- filter(immdata_imms, Indiv %in% c(pedigree$Dad, pedigree$Mom))

# How many immigrants are in each cohort?
num_MF_imms_per_cohort <- left_join(immbreeders, pedigree) %>%
  group_by(ImmCohort, Sex) %>%
  dplyr::summarize(num_imms = n()) %>% 
  filter(ImmCohort < 1998)

####Part D: Autosomal exp gen contribs 15 years after immigration vs imm cohort size####

#get data
obsImmYrA<-read.table('working_files/intermediate_files/ImmContribYearly_malevsfemale_A.1.drop.data.txt',header=TRUE)
simImmYrA<-read.table('working_files/intermediate_files/ImmContribYearly_malevsfemale_A.1.drop.sim.txt',header=TRUE)

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

# convert allele to cohort year
# males have allele 2 = 1991, allele 4 = 1992, allele 6 = 1993, etc
# females have allele 3 = 1991, allele 5 = 1992, allele 7 = 1993, etc
# so for males it is allele/2 + 1990
# and for females it is (allele-1)/2 + 1990
simsumImmYrAFplus15 <- mutate(simsumImmYrAF, ImmCohort = ((allele-1)/2) + 1990, Sex = 2) %>%
  filter(Year==ImmCohort+15, ImmCohort < 1998)
simsumImmYrAMplus15 <- mutate(simsumImmYrAM, ImmCohort = (allele/2) + 1990, Sex = 1) %>%
  filter(Year==ImmCohort+15, ImmCohort < 1998)

simsumImmYrAplus15 <- bind_rows(simsumImmYrAFplus15, simsumImmYrAMplus15) %>%
  left_join(num_MF_imms_per_cohort, by = c("ImmCohort", "Sex")) %>%
  mutate(Sex = ifelse(Sex == 1, "M", "F"))

# plot
(A_imm_contribs_vs_cohort_size <- ggplot(simsumImmYrAplus15, aes(x=num_imms, y=mean, color=Sex)) +
  geom_point(show.legend = FALSE) +
  geom_linerange(aes(ymin=q1, ymax=q3), show.legend = FALSE) +
  geom_smooth(method = "lm", aes(fill = Sex), alpha = 0.25, size = 0.5) +
  scale_color_manual(values = c("indianred1", "cornflowerblue")) +
  scale_fill_manual(values = c("indianred1", "cornflowerblue")) +
  labs(x="Size of immigrant cohort, sex-specific \n(number of birds of the specified sex)", y="Autosomal expected genetic \ncontrib. after 15 years") +
  plottheme +
  guides(color=guide_legend(override.aes=list(color=NA))) +
  theme(legend.position = "none",plot.margin=unit(c(0.2,0.1,0,0.15),'cm')))

simsumImmYrAplus15_lm <- lm(data=simsumImmYrAplus15, mean ~ num_imms + Sex)
summary(simsumImmYrAplus15_lm)

####Part E: Z exp gen contribs 15 years after immigration vs imm cohort size####

#get data
obsImmYrZ<-read.table('working_files/intermediate_files/ImmContribYearly_malevsfemale_Z.1.drop.data.txt',header=TRUE)
simImmYrZ<-read.table('working_files/intermediate_files/ImmContribYearly_malevsfemale_Z.1.drop.sim.txt',header=TRUE)

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
simImmYrZM$totAlleles<-lapply(c(1:length(simImmYrZM$cohort)), function(x) unique(obsImmYrZ[obsImmYrZ$cohort_year==simImmYrZM[x,'Year'],'allele_count']))

#get mean 
simImmYrZAvgM<-ddply(simImmYrZM,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count,all_alleles_count))))

simImmYrZM$allele_count <- as.numeric(simImmYrZM$allele_count)
simImmYrZM$totAlleles <- as.numeric(simImmYrZM$totAlleles)


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
simImmYrZF$totAlleles<-lapply(c(1:length(simImmYrZF$cohort)), function(x) unique(obsImmYrZ[obsImmYrZ$cohort_year==simImmYrZF[x,'Year'],'allele_count']))

#get mean 
simImmYrZAvgF<-ddply(simImmYrZF,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count,all_alleles_count))))

simImmYrZF$allele_count <- as.numeric(simImmYrZF$allele_count)
simImmYrZF$totAlleles <- as.numeric(simImmYrZF$totAlleles)

simsumImmYrZF<-ddply(simImmYrZF,.(Year,allele),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))
simsumImmYrZF<-simsumImmYrZF[simsumImmYrZF$allele != 1,]

# convert allele to cohort year
# males have allele 2 = 1991, allele 4 = 1992, allele 6 = 1993, etc
# females have allele 3 = 1991, allele 5 = 1992, allele 7 = 1993, etc
# so for males it is allele/2 + 1990
# and for females it is (allele-1)/2 + 1990
simsumImmYrZFplus15 <- mutate(simsumImmYrZF, ImmCohort = ((allele-1)/2) + 1990, Sex = 2) %>%
  filter(Year==ImmCohort+15, ImmCohort < 1998)
simsumImmYrZMplus15 <- mutate(simsumImmYrZM, ImmCohort = (allele/2) + 1990, Sex = 1) %>%
  filter(Year==ImmCohort+15, ImmCohort < 1998)

simsumImmYrZplus15 <- bind_rows(simsumImmYrZFplus15, simsumImmYrZMplus15) %>%
  left_join(num_MF_imms_per_cohort, by = c("ImmCohort", "Sex")) %>%
  mutate(Sex = ifelse(Sex == 1, "M", "F"))

# plot
(Z_imm_contribs_vs_cohort_size <- ggplot(simsumImmYrZplus15, aes(x=num_imms, y=mean, color=Sex)) +
  geom_point(show.legend = FALSE) +
  geom_linerange(aes(ymin=q1, ymax=q3), show.legend = FALSE) +
  geom_smooth(method = "lm", aes(fill = Sex), alpha = 0.25, size = 0.5) +
  scale_color_manual(values = c("indianred1", "cornflowerblue")) +
  scale_fill_manual(values = c("indianred1", "cornflowerblue")) +
  labs(x="Size of immigrant cohort, sex-specific \n(number of birds of the specified sex)", y="Z expected genetic \ncontrib. after 15 years") +
  plottheme +
  guides(color=guide_legend(override.aes=list(color=NA))) +
  theme(legend.position = "none",plot.margin=unit(c(0.2,0.1,0,0.15),'cm')))

simsumImmYrAplus15_lm <- lm(data=simsumImmYrZplus15, mean ~ num_imms + Sex)
summary(simsumImmYrAplus15_lm)


####Putting it all together####

(left_side <- plot_grid(sexed_breed_nonbreed, unsexed_breed_nonbreed, A_imm_contribs_vs_cohort_size, ncol = 1, nrow = 3, rel_heights = c(0.75, 0.25, 0.5), align = "hv", axis = "lrtb", labels = c("A", "", "D")))
dev.off()
(right_side <- plot_grid(A_imm_gen_contribs, Z_imm_gen_contribs, Z_imm_contribs_vs_cohort_size, ncol = 1, nrow = 3, rel_heights = c(0.5, 0.5, 0.5), align = "hv", axis = "lrtb", labels = c("B", "C", "E")))
dev.off()
(fig3_no_legend <- plot_grid(left_side, right_side, ncol = 2, nrow = 1, rel_widths = c(0.5,0.5), align = "hv", axis = "lrtb"))
plot_grid(sexed_unsexed_breed_nonbreed_legend, fig3_no_legend, ncol = 1, nrow = 2, rel_heights = c(0.05,1))
#ggsave("fig3_immcontrib.pdf", width = 6.5, height = 6.5, units = "in")


####Supplemental figure 4: Z/A for male & female imms####

# combine A and Z data into a single table for plotting
simsumImm_from1990 <- select(simsumImmA_from1990, year, allele, mean_A = mean, q1_A = q1, q3_A = q3) %>%
  left_join(select(simsumImmZ_from1990, year, allele, mean_Z = mean, q1_Z = q1, q3_Z = q3))

ggplot(simsumImm_from1990,aes(y=mean_Z/mean_A,x=year)) +
  annotate('segment', x=1993,xend=2013,y=1,yend=1,alpha=0.1)+
  geom_hline(yintercept = 4/3, linetype="dashed",alpha=0.5)+
  annotate('text',x=1991,y=1.4,label='4/3',size=3)+
  geom_hline(yintercept = 2/3, linetype="dashed",alpha=0.5)+
  annotate('text',x=1991,y=0.6,label='2/3',size=3)+
  geom_line(aes(color=as.factor(allele)), show.legend = FALSE)+
  geom_ribbon(aes(x = year, ymin = q1_Z/mean_A, ymax = q3_Z/mean_A, fill = as.factor(allele)), alpha=0.3) +
  geom_ribbon(aes(x = year, ymin = q1_A/mean_A, ymax = q3_A/mean_A, fill = as.factor(allele)), alpha=0.3) +
  scale_color_manual(labels = c("Male", "Female"), values = c("cornflowerblue", "indianred1")) +
  scale_fill_manual(labels = c("Male", "Female"), values = c("cornflowerblue", "indianred1")) +
  labs(x = "Year", y = "Z/A Immigrant Contributions",color="Sex",fill="Sex") +
  plottheme + 
  theme(legend.key.size = unit(0.2, "in"))
#ggsave("ZA_ratio_immigrants.pdf", width = 5, height = 3.5, units = "in")  

