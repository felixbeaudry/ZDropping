#script to make input files for alleleFreqModel from plink and pedigree
#Felix Beaudry & Rose Driscoll 
#Last updated: Oct 15 2021

library(tidyverse)
library(data.table)
'%ni%' <- Negate('%in%')

####top####

ped_geno <- fread('genotypeFiles/geno.anon.ped')
ped <- ped_geno[,c(1:6)]

write.table(ped, 
            file = "alleleFreqModel/working_files/pedigree.txt", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, qmethod = c("escape", "double"))

map<-read.table('genotypeFiles/geno.map')
map$map_pos <- seq(1,length(map$V1))
write.table(map, 
            file = "alleleFreqModel/working_files/geno.map", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, qmethod = c("escape", "double"))

#calculate positions for data subsets
#Apos <- sort(c(5+(2*map$map_pos[map$V1 %in% c(0:38)]),6+(2*map$map_pos[map$V1 %in% c(0:38)])))
Zpos <- sort(c(5+(2*map$map_pos[map$V1 == 39]),6+(2*map$map_pos[map$V1 == 39])))
PARpos <- sort(c(5+(2*map$map_pos[map$V1 == 41]),6+(2*map$map_pos[map$V1 == 41])))

#make positions subsets - for selection gene dropping
zpedgeno <- cbind.data.frame(ped_geno[,c(1:6)],ped_geno[,..Zpos])
PARpedgeno <- cbind.data.frame(ped_geno[,c(1:6)],ped_geno[,..PARpos])

write.table(zpedgeno, 
            file = "GeneDropping/working_files/intermediate_files/FSJpedgeno_Zsexlinked.ped", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))

write.table(PARpedgeno, 
            file = "GeneDropping/working_files/intermediate_files/FSJpedgeno_Zpseudoautosomal.ped", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))

####make pedgenos for allelefreq####
indivlist <- fread('genotypeFiles/indivlist.txt')

#add assigned sexes of unsexed birds back to indivlist 
#(this way, a given unsexed bird will always have the same assigned sex even if it appears multiple times in indivlist)

#check for unsexed indivs & assign them a sex
#unsexed_indivs <- ped$Indiv[ped$Sex==0]
#simulated_sexes <- sample(x = c(1,2), size = length(unsexed_indivs), prob = c(0.5,0.5), replace = TRUE)
#ped$simSex <- ped$Sex
#ped$simSex[ped$simSex==0] <- simulated_sexes

#indivlist <- left_join(indivlist,ped,by=c("USFWS"="Indiv"))
#names(indivlist)<-c('Year','Indiv','Category','Genotyped','Mom','Dad','Sex','simSex')

indivlist$Indiv <- as.character(indivlist$Indiv)

#
ped_genod <- ped_geno[ped_geno$ID %in% indivlist$Indiv[indivlist$Genotyped == "Y"],]

tmpped<-ped_genod
ped_genod<-tmpped[,2:5]

#convert genotype data to one column per SNP
for (x in seq(7,(6 + length(map$V2)*2),by=2)) { #cycle through SNPs
  cols <- c(x,x+1)
  thisped<-tmpped[,..cols]
  thisped$geno<-rowSums(thisped[,1:2])-2
  thisped[thisped$geno<0,'geno']<-NA
  ped_genod<-cbind(ped_genod,thisped$geno)
}

names(ped_genod)<-c('Indiv','Dad','Mom','Sex',paste("SNP",c(1:(length(map$V2))),sep=""))

#correct female genotypes on Z to indicate that they are heterozygous (1 instead of 2)
Zpos_g <- c(4 + map$map_pos[map$V1 == 39])
ped_genod[ped_genod$Sex==2,Zpos_g] <- ped_genod[ped_genod$Sex==2,..Zpos_g]/2

ped_genod$Indiv <- as.character(ped_genod$Indiv)

genocols <- c(1,5:length(ped_genod))
indivlistgeno <- left_join(indivlist,ped_genod[,..genocols],by=c("Indiv"="Indiv"))

Apos_g <- c(8 + map$map_pos[map$V1 %in% c(0:38)])
Zpos_g <- c(8 + map$map_pos[map$V1 == 39])

indivlistgeno_A <- cbind.data.frame(indivlistgeno[,c(1:8)],indivlistgeno[,..Apos_g])
indivlistgeno_Z <- cbind.data.frame(indivlistgeno[,c(1:8)],indivlistgeno[,..Zpos_g])

save(indivlistgeno_A,file='alleleFreqModel/working_files/intermediate_files/indivlistgeno_A.rdata')
save(indivlistgeno_Z,file='alleleFreqModel/working_files/intermediate_files/indivlistgeno_Z.rdata')
