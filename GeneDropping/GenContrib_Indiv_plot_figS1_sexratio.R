#29 March 2024
#Rose Driscoll
#script to plot the genetic contributions of male and female breeders for autosomes and Z with unsexed birds assigned a sex with sex ratio 0, 0.5, or 1
#just the plotting commands

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
breeders_926$Indiv <- as.character(breeders_926$Indiv)

# Read in pedigree
pedigree <- read.table("working_files/pedigree.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
colnames(pedigree) <- c("Fam", "Indiv", "Dad", "Mom", "Sex", "Pheno")
pedigree$Indiv <- as.character(pedigree$Indiv)
pedigree$Mom <- as.character(pedigree$Mom)
pedigree$Dad <- as.character(pedigree$Dad)


## Get 2013 Z genetic contributions (sex ratio 0, 0.5, 1) for each individual

# Make a table to save individuals' 2013 contributions in
indiv_contribs_A_Z_sexratio_2013<-NULL

# Only running this for breeders that actually have kids
breeders_926$Kids <- breeders_926$Indiv %in% c(pedigree$Dad, pedigree$Mom)

# Go through all of the individuals
for (i in breeders_926[breeders_926$Kids,'Indiv']) {
  
  #get autosomal, sex ratio 0 data for individual i
  obsIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.A.0.1.drop.data.txt',sep=''),header=TRUE)
  simIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.A.0.1.drop.sim.txt',sep=''),header=TRUE)
  
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
  
  #save mean autosomal contribution (sex ratio 0) for 2013
  indivauto0 <- mutate(simIndivSum, Indiv = i) %>%
    select(Indiv, auto_0_mean = mean)
  
  
  #get autosomal, sex ratio 0.5 data for individual i
  obsIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.A.0.5.1.drop.data.txt',sep=''),header=TRUE)
  simIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.A.0.5.1.drop.sim.txt',sep=''),header=TRUE)
  
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
  
  #save mean autosomal contribution (sex ratio 0.5) for 2013
  indivauto0.5 <- mutate(simIndivSum, Indiv = i) %>%
    select(Indiv, auto_0.5_mean = mean) 
  
  #get autosomal, sex ratio 1 data for individual i
  obsIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.A.1.1.drop.data.txt',sep=''),header=TRUE)
  simIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.A.1.1.drop.sim.txt',sep=''),header=TRUE)
  
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
  
  #save mean autosomal contribution (sex ratio 1) for 2013
  indivauto1 <- mutate(simIndivSum, Indiv = i) %>%
    select(Indiv, auto_1_mean = mean)
  
  
  #get Z, sex ratio 0 data for individual i
  obsIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.Z.0.1.drop.data.txt',sep=''),header=TRUE)
  simIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.Z.0.1.drop.sim.txt',sep=''),header=TRUE)
  
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
  
  #save mean Z contribution (sex ratio 0) for 2013
  indivZ0 <- mutate(simIndivSum, Indiv = i) %>%
    select(Indiv, Z_0_mean = mean)
  
  
  #get Z, sex ratio 0.5 data for individual i
  obsIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.Z.0.5.1.drop.data.txt',sep=''),header=TRUE)
  simIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.Z.0.5.1.drop.sim.txt',sep=''),header=TRUE)
  
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
  
  #save mean Z contribution (sex ratio 0.5) for 2013
  indivZ0.5 <- mutate(simIndivSum, Indiv = i) %>%
    select(Indiv, Z_0.5_mean = mean)
  
  
  #get Z, sex ratio 1 data for individual i
  obsIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.Z.1.1.drop.data.txt',sep=''),header=TRUE)
  simIndiv<-read.table(file=paste('working_files/intermediate_files/IndivContrib_',i,'.ped.Z.1.1.drop.sim.txt',sep=''),header=TRUE)
  
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
  
  #save mean Z contribution (sex ratio 1) for 2013
  indivZ1 <- mutate(simIndivSum, Indiv = i) %>%
    select(Indiv, Z_1_mean = mean)
  
  
  # add individual to table
  indivrow <- right_join(indivauto0, indivauto0.5) %>%
    right_join(indivauto1) %>%
    right_join(indivZ0) %>%
    right_join(indivZ0.5) %>%
    right_join(indivZ1)
  indiv_contribs_A_Z_sexratio_2013<-rbind(indiv_contribs_A_Z_sexratio_2013, indivrow)
  
}
# Add back in breeders with no kids, with 0 contributions
indiv_contribs_A_Z_sexratio_2013 <- left_join(select(breeders_926, Indiv), indiv_contribs_A_Z_sexratio_2013,by="Indiv") %>%
  mutate(auto_0_mean = ifelse(is.na(auto_0_mean), 0, auto_0_mean), 
         auto_0.5_mean = ifelse(is.na(auto_0.5_mean), 0, auto_0.5_mean), 
         auto_1_mean = ifelse(is.na(auto_1_mean), 0, auto_1_mean),
         Z_0_mean = ifelse(is.na(Z_0_mean), 0, Z_0_mean), 
         Z_0.5_mean = ifelse(is.na(Z_0.5_mean), 0, Z_0.5_mean), 
         Z_1_mean = ifelse(is.na(Z_1_mean), 0, Z_1_mean))

# add sex info from pedigree
indiv_contribs_A_Z_sexratio_2013_sex <- left_join(indiv_contribs_A_Z_sexratio_2013, select(pedigree, Indiv, Sex))
indiv_contribs_A_Z_sexratio_2013_sex$Sex <- as.factor(indiv_contribs_A_Z_sexratio_2013_sex$Sex)


## Plot

# Plot comparing 0 vs 1
ggplot(indiv_contribs_A_Z_sexratio_2013_sex, aes(x=auto_0_mean, y=auto_1_mean, color=Sex)) +
  geom_point(size = 1.2) +
  scale_color_manual(values = c("cornflowerblue", "indianred1")) +
  plottheme +
  labs(x="Autosomal expected genetic contrib. in 2013\nSex ratio 0 applied to unsexed nestlings", y="Autosomal expected genetic contrib. in 2013\nSex ratio 1 applied to unsexed nestlings") +
  theme(legend.position = "none") ->auto_0_vs_1
ggplot(indiv_contribs_A_Z_sexratio_2013_sex, aes(x=Z_0_mean, y=Z_1_mean, color=Sex)) +
  geom_point(size = 1.2) +
  scale_color_manual(values = c("cornflowerblue", "indianred1")) +
  plottheme +
  labs(x="Z expected genetic contrib. in 2013\nSex ratio 0 applied to unsexed nestlings", y="Z expected genetic contrib. in 2013\nSex ratio 1 applied to unsexed nestlings") +
  theme(legend.position = "none")->Z_0_vs_1
pdf("fig_S1_EGC_AZ_sexratio0v1.pdf",width=6.5,height=3.25)
plot_grid(auto_0_vs_1, Z_0_vs_1,labels = "AUTO")
dev.off()

# Correlation test for the above
cor.test(indiv_contribs_A_Z_sexratio_2013_sex$auto_0_mean, 
         indiv_contribs_A_Z_sexratio_2013_sex$auto_1_mean, method = "spearman")
cor.test(indiv_contribs_A_Z_sexratio_2013_sex$Z_0_mean, 
         indiv_contribs_A_Z_sexratio_2013_sex$Z_1_mean, method = "spearman")

