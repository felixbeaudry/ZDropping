#10 March 2021
#Rose Driscoll
#script to plot the pedigree, genetic & genealogical contributions of a male and female breeder 
#for autosomes and Z


## Setup

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(kinship2)
library(cowplot)
library(gridGraphics)

# Plot theme
plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3), 
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"), 
                    plot.background = element_rect(fill = "white"),  
                    axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), 
                    axis.title = element_text(size=7), plot.title = element_text(size=8), 
                    legend.position="right", legend.text = element_text(size=7),
                    legend.title = element_text(size=8), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(1,"cm"))

## Read in data

# Read in 1990-2013 breeders data
load("indivGenContrib.rdata")

# Read in and tidy pedigree
ped <- read.table("FSJpedgeno_Zsexlinked.ped", header = FALSE, sep = " ", stringsAsFactors = FALSE)
pedigree <- ped[,1:6]
colnames(pedigree) <- c("Fam", "Indiv", "Dad", "Mom", "Sex", "Pheno")

# Read in all of the data tables including breeders data
load("all_tables_v2.rdata")

## Genealogical contributions for all breeders

# Based on the pedigree, make a list of the descendants of every breeder we'll be looking at
descendants <- NULL
for (i in dead2013age1$USFWS) {
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

# Read in data on individuals including natal years
indiv_data<-read.table('IndivDataUSFWS.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
nestling_data<-select(indiv_data, Indiv, NatalYear)

# add natal year info to descendants dataset
descendants_years <- left_join(descendants, nestling_data, by = c("descendant_id" = "Indiv"))

# how many nestlings are there each year?
num_nestlings_yearly <- group_by(nestling_data, NatalYear) %>%
  dplyr::summarize(num_nestlings = n())


## Figure 1: genetic (A & Z) and genealogical contributions across years for a mated pair

# The pair is:
m<-"1513-50363"
f<-"1043-90505"

# make a table for this pair's genealogical contributions by year
# fill out the years where they don't have nestlings with 0s
filter(descendants_years, indiv_id == "1513-50363") %>%
  group_by(NatalYear) %>%
  dplyr::summarize(num_descendants = n()) %>%
  right_join(num_nestlings_yearly) %>%
  mutate(num_descendants = ifelse(is.na(num_descendants), 0, num_descendants), 
         prop_nestlings = num_descendants/num_nestlings) -> genealogical_contribs_by_year


# Part A: the pedigree

# Get a list of all the individuals we want in this pedigree
# they have all the same descendants so we can just get one's descendants
focal_family <- filter(descendants, indiv_id == f)
# add the two focal individuals to the list
focal_family_indivs <- c(m, f, focal_family$descendant_id)

# Get the pedigree data for this family and modify it a bit for the plot
focal_family_insiders <- filter(pedigree, Indiv %in% focal_family_indivs) %>%
  select(Indiv, Dad, Mom, Sex) %>%
  mutate(Mom = ifelse(Indiv==m | Indiv==f, NA, Mom), Dad = ifelse(Indiv==m | Indiv == f, NA, Dad))
# Also need the second parent for individuals who have only 1 parent in this tree (make these founders)
focal_family_outsiders <- filter(pedigree, Indiv %in% focal_family_insiders$Mom | 
                                   Indiv %in% focal_family_insiders$Dad, !(Indiv %in% focal_family_indivs)) %>%
  mutate(Mom = NA, Dad = NA)
# Add those extra parents/outsiders to the family pedigree
focal_family_ped <- bind_rows(focal_family_insiders, focal_family_outsiders) %>%
  mutate(Sex = ifelse(Sex==0, 3, Sex)) %>%
  select(Indiv, Dad, Mom, Sex)
# We will also need a status variable that says whether the indiv is alive in 2013
# Get last observation for each indiv
focal_family_ped_status <- filter(all.obs, USFWS %in% focal_family_ped$Indiv) %>%
  group_by(USFWS) %>%
  dplyr::summarize(last_year = max(Year)) %>%
  mutate(Status = ifelse(last_year > 2013, 1, 0)) %>%
  select(Indiv = USFWS, Status) %>%
  right_join(focal_family_ped)

# Draw the pedigree w/ kinship2 package
focal_family_ped_to_plot <- with(focal_family_ped_status[,c(1,3:5,2)], pedigree(Indiv, Dad, Mom, Sex, Status))
plot(focal_family_ped_to_plot, id=rep("", times = nrow(focal_family_ped)), symbolsize = 1)
# use recordPlot() to capture the base R plot in a variable
mated_pair_ped <- recordPlot()


# Part B: autosomal expected genetic contributions vs. genealogical contributions

#get autosomal expected genetic contribution data for the male
obsIndivMA<-read.table(file=paste('IndivContrib_',m,'.ped.A.1.drop.data.txt',sep=''),header=TRUE)
simIndivMA<-read.table(file=paste('IndivContrib_',m,'.ped.A.1.drop.sim.txt',sep=''),header=TRUE)

#make sure number of simulations adds up to 1000000
ddply(simIndivMA[simIndivMA$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))
for (y in c(0:23)){
  if (sum(simIndivMA[simIndivMA$allele==2 & simIndivMA$cohort==y,'all_alleles_count'])<1000000) {
    simIndivMA<-rbind(simIndivMA,cbind(cohort=y,allele=2,allele_count=0,all_alleles_count=1000000-sum(simIndivMA[simIndivMA$allele==2 & simIndivMA$cohort==y,'all_alleles_count'])))
  }
}
ddply(simIndivMA[simIndivMA$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))

#add year info
simIndivMA$Year<-simIndivMA$cohort + 1990

#subset data to only include 1990-2013 and allele 2 (i.e., alleles marked as coming from this indiv)
simIndivMA<-simIndivMA[simIndivMA$allele==2 & simIndivMA$Year < 2014,]

#get total number of alleles sampled each year
simIndivMA$totAlleles<-laply(c(1:nrow(simIndivMA)), function(x) unique(obsIndivMA[obsIndivMA$cohort_year==simIndivMA[x,'Year'],'all_alleles_count']))

#summarize data to get the mean, q1, and q3 of the simulations
simIndivMASum<-ddply(simIndivMA,.(Year),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), 
                   q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), 
                   q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))


#get autosomal expected genetic contribution data for the female
obsIndivFA<-read.table(file=paste('IndivContrib_',f,'.ped.A.1.drop.data.txt',sep=''),header=TRUE)
simIndivFA<-read.table(file=paste('IndivContrib_',f,'.ped.A.1.drop.sim.txt',sep=''),header=TRUE)

#make sure number of simulations adds up to 1000000
ddply(simIndivFA[simIndivFA$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))
for (y in c(0:23)){
  if (sum(simIndivFA[simIndivFA$allele==2 & simIndivFA$cohort==y,'all_alleles_count'])<1000000) {
    simIndivFA<-rbind(simIndivFA,cbind(cohort=y,allele=2,allele_count=0,all_alleles_count=1000000-sum(simIndivFA[simIndivFA$allele==2 & simIndivFA$cohort==y,'all_alleles_count'])))
  }
}
ddply(simIndivFA[simIndivFA$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))

#add year info
simIndivFA$Year<-simIndivFA$cohort + 1990

#subset data to only include 1990-2013 and allele 2 (i.e., alleles marked as coming from this indiv)
simIndivFA<-simIndivFA[simIndivFA$allele==2 & simIndivFA$Year < 2014,]

#get total number of alleles sampled each year
simIndivFA$totAlleles<-laply(c(1:nrow(simIndivFA)), function(x) unique(obsIndivFA[obsIndivFA$cohort_year==simIndivFA[x,'Year'],'allele_count']))

#summarize data to get the mean, q1, and q3 of the simulations
simIndivFASum<-ddply(simIndivFA,.(Year),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), 
                   q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), 
                   q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))


# We only want contributions to descendants, but the male is counted the year he is born (1999)
# Delete these counts as they are not representative of a descendant.
simIndivMASum[simIndivMASum$Year==1999,"mean"]<-0
simIndivMASum[simIndivMASum$Year==1999,"q1"]<-0
simIndivMASum[simIndivMASum$Year==1999,"q3"]<-0


# plot autosomal expected genetic contributions + genealogical contributions
auto_plot <- ggplot() + plottheme + xlab('Year') + ylab('Genealogical & \nautosomal expected genetic contribution') +
  geom_ribbon(data=simIndivFASum,aes(x=Year,ymin=q1,ymax=q3),fill='mediumpurple', alpha = 0.15) +
  geom_line(data=simIndivFASum,aes(x=Year, y=mean),size=0.5,col='mediumpurple') +
  geom_ribbon(data = simIndivMASum,aes(x=Year,ymin=q1,ymax=q3),fill='mediumpurple',alpha=0.15) +
  geom_line(data=simIndivMASum,aes(x=Year, y=mean),size=0.5,col='mediumpurple') +
  geom_line(data = genealogical_contribs_by_year, aes(x=NatalYear, y=prop_nestlings),size=0.5,col='gray50', linetype = "dashed")


# Part C: Z expected genetic contributions vs. genealogical contributions

#get Z expected genetic contribution data for the male

obsIndivMZ<-read.table(file=paste('IndivContrib_',m,'.ped.Z.1.drop.data.txt',sep=''),header=TRUE)
simIndivMZ<-read.table(file=paste('IndivContrib_',m,'.ped.Z.1.drop.sim.txt',sep=''),header=TRUE)

#make sure number of simulations adds up to 1000000
ddply(simIndivMZ[simIndivMZ$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))
for (y in c(0:23)){
  if (sum(simIndivMZ[simIndivMZ$allele==2 & simIndivMZ$cohort==y,'all_alleles_count'])<1000000) {
    simIndivMZ<-rbind(simIndivMZ,cbind(cohort=y,allele=2,allele_count=0,all_alleles_count=1000000-sum(simIndivMZ[simIndivMZ$allele==2 & simIndivMZ$cohort==y,'all_alleles_count'])))
  }
}
ddply(simIndivMZ[simIndivMZ$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))

#add year info
simIndivMZ$Year<-simIndivMZ$cohort + 1990

#subset data to only include 1990-2013 and allele 2 (i.e., alleles marked as coming from this indiv)
simIndivMZ<-simIndivMZ[simIndivMZ$allele==2 & simIndivMZ$Year < 2014,]

#get total number of alleles sampled each year
simIndivMZ$totAlleles<-laply(c(1:nrow(simIndivMZ)), function(x) unique(obsIndivMZ[obsIndivMZ$cohort_year==simIndivMZ[x,'Year'],'all_alleles_count']))

#summarize data to get the mean, q1, and q3 of the simulations
simIndivMZSum<-ddply(simIndivMZ,.(Year),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), 
                   q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), 
                   q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))

#get Z expected genetic contribution data for the female
obsIndivFZ<-read.table(file=paste('IndivContrib_',f,'.ped.Z.1.drop.data.txt',sep=''),header=TRUE)
simIndivFZ<-read.table(file=paste('IndivContrib_',f,'.ped.Z.1.drop.sim.txt',sep=''),header=TRUE)

#make sure number of simulations adds up to 1000000
ddply(simIndivFZ[simIndivFZ$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))
for (y in c(0:23)){
  if (sum(simIndivFZ[simIndivFZ$allele==2 & simIndivFZ$cohort==y,'all_alleles_count'])<1000000) {
    simIndivFZ<-rbind(simIndivFZ,cbind(cohort=y,allele=2,allele_count=0,all_alleles_count=1000000-sum(simIndivFZ[simIndivFZ$allele==2 & simIndivFZ$cohort==y,'all_alleles_count'])))
  }
}
ddply(simIndivFZ[simIndivFZ$allele==2,],.(cohort),summarize, sum=sum(all_alleles_count))

#add year info
simIndivFZ$Year<-simIndivFZ$cohort + 1990

#subset data to only include 1990-2013 and allele 2 (i.e., alleles marked as coming from this indiv)
simIndivFZ<-simIndivFZ[simIndivFZ$allele==2 & simIndivFZ$Year < 2014,]

#get total number of alleles sampled each year
simIndivFZ$totAlleles<-laply(c(1:nrow(simIndivFZ)), function(x) unique(obsIndivFZ[obsIndivFZ$cohort_year==simIndivFZ[x,'Year'],'allele_count']))

#summarize data to get the mean, q1, and q3 of the simulations
simIndivFZSum<-ddply(simIndivFZ,.(Year),summarize, mean=mean(unlist(rep(allele_count/totAlleles,all_alleles_count))), 
                   q1=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.025), 
                   q3=quantile(unlist(rep(allele_count/totAlleles,all_alleles_count)),0.975))

# We only want contributions to descendants, but the male is counted the year he is born (1999)
# Delete these counts as they are not representative of a descendant.
simIndivMZSum[simIndivMASum$Year==1999,"mean"]<-0
simIndivMZSum[simIndivMASum$Year==1999,"q1"]<-0
simIndivMZSum[simIndivMASum$Year==1999,"q3"]<-0


# plot Z expected genetic contributions + genealogical contributions
Z_plot <- ggplot() + plottheme + xlab('Year') + ylab('Genealogical & \nZ genetic contribution') +
  geom_ribbon(data=simIndivFZSum,aes(x=Year,ymin=q1,ymax=q3),fill='indianred1', alpha = 0.3) +
  geom_line(data=simIndivFZSum,aes(x=Year, y=mean),size=0.5,col='indianred1') +
  geom_ribbon(data=simIndivMZSum,aes(x=Year,ymin=q1,ymax=q3),fill='cornflowerblue', alpha = 0.3) +
  geom_line(data=simIndivMZSum,aes(x=Year, y=mean),size=0.5,col='cornflowerblue') +
  geom_line(data = genealogical_contribs_by_year, aes(x=NatalYear, y=prop_nestlings),size=0.5,col='gray50', linetype = "dashed")


## Combine plots

# add some space at the top of these 2 plots to make room for the letters
auto_margin <- auto_plot + theme(plot.margin = unit(c(10,5.5,5.5,5.5), "pt"))
Z_margin <- Z_plot + theme(plot.margin = unit(c(10,5.5,5.5,5.5), "pt"))

plot_grid(mated_pair_ped, auto_margin, Z_margin, nrow = 1, ncol = 3, labels = "AUTO")
#ggsave("plotIndivGenContrib_tidy_20210310/fig1_ped_A_Z.pdf", width = 6.5, height = 2, units = "in")

