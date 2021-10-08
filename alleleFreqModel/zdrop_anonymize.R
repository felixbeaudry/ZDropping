
library(tidyverse)
'%ni%' <- Negate('%in%')

####top####

ped_only <- fread('working_files/pedigree.txt')
maleCode <- unique(ped_only$Sex[ped_only$Indiv %in% ped_only$Dad]) #makes var that contains what number is male
femaleCode <- unique(ped_only$Sex[ped_only$Indiv %in% ped_only$Mom]) #makes var that contains what number is female

write.table(ped_only[,c(1,2)], 
            file = "working_files/intermediate_files/keepInds.list", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))

###
#run plink 
###

ped<-read.table('working_files/FSJfullPedFiltDogFINAL12July2016finalSexNumMAF05geno.ped',header=FALSE,sep=' ',stringsAsFactors=FALSE)
ped_long <- left_join(ped_only,ped[,-c(1,3:6)],by=c("Indiv"="V2")) #add in "placeholder individuals that we added during our demographic data QA/QC process (after the genotype data were finalized)"
ped_long[is.na(ped_long)] <- 0

map<-read.table('working_files/FSJfullPedFiltDogFINAL12July2016finalSexNumMAF05geno.map')
map$map_pos <- seq(1,length(map$V1))

#calculate positions for data subsets
#Apos <- sort(c(5+(2*map$map_pos[map$V1 %in% c(0:38)]),6+(2*map$map_pos[map$V1 %in% c(0:38)])))
Zpos <- sort(c(5+(2*map$map_pos[map$V1 == 39]),6+(2*map$map_pos[map$V1 == 39])))
PARpos <- sort(c(5+(2*map$map_pos[map$V1 == 41]),6+(2*map$map_pos[map$V1 == 41])))

#make positions subsets - for selection gene dropping
zpedgeno <- cbind.data.frame(ped_long[,c(1:6)],ped_long[,..Zpos])
PARpedgeno <- cbind.data.frame(ped_long[,c(1:6)],ped_long[,..PARpos])

write.table(zpedgeno, 
            file = "working_files/intermediate_files/FSJpedgeno_Zsexlinked.ped", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))

write.table(PARpedgeno, 
            file = "working_files/intermediate_files/FSJpedgeno_Zpseudoautosomal.ped", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))

####make pedgenos for allelefreq####
load('working_files/simindivFIXmin2obs.rdata')

#add assigned sexes of unsexed birds back to indivlist 
#(this way, a given unsexed bird will always have the same assigned sex even if it appears multiple times in indivlist)

#check for unsexed indivs & assign them a sex
#unsexed_indivs <- ped_only$Indiv[ped_only$Sex==0]
#simulated_sexes <- sample(x = c(1,2), size = length(unsexed_indivs), prob = c(0.5,0.5), replace = TRUE)
#ped_only$simSex <- ped_only$Sex
#ped_only$simSex[ped_only$simSex==0] <- simulated_sexes

simsex <- fread('working_files/FSJ_sex_data_real_and_simulated_20201015.csv')
names(simsex) <- c("USFWS","simSex")

indivlist <- left_join(simindivFIXmin2obs,ped_only[,c(2,5)],by=c("USFWS"="Indiv")) %>% left_join(simsex)



names(indivlist)<-c('Year','Indiv','Category','Genotyped','Mom','Dad','Sex','simSex')

indivlist$Indiv <- as.character(indivlist$Indiv)

#
ped_genod <- ped_long[ped_long$Indiv %in% simindivFIXmin2obs$USFWS[simindivFIXmin2obs$genotyped == "Y"],]

tmpped<-ped_genod
pedgeno<-tmpped[,2:5]

#convert genotype data to one column per SNP
for (x in seq(7,(6 + length(map$V2)*2),by=2)) { #cycle through SNPs
  cols <- c(x,x+1)
  thisped<-tmpped[,..cols]
  thisped$geno<-rowSums(thisped[,1:2])-2
  thisped[thisped$geno<0,'geno']<-NA
  pedgeno<-cbind(pedgeno,thisped$geno)
}

names(pedgeno)<-c('Indiv','Dad','Mom','Sex',paste("SNP",c(1:(length(map$V2))),sep=""))

#correct female genotypes on Z to indicate that they are heterozygous (1 instead of 2)
Zpos_g <- c(4 + map$map_pos[map$V1 == 39])
pedgeno[pedgeno$Sex==2,Zpos_g] <- pedgeno[pedgeno$Sex==2,..Zpos_g]/2

#save(pedgeno,file="working_files/intermediate_files/pedgeno.rdata")

pedgeno$Indiv <- as.character(pedgeno$Indiv)

genocols <- c(1,5:length(pedgeno))
indivlistgeno <- left_join(indivlist,pedgeno[,..genocols],by=c("Indiv"="Indiv"))

#save(indivlistgeno,file='working_files/intermediate_files/indivlist_geno.rdata')

Apos_g <- c(8 + map$map_pos[map$V1 %in% c(0:38)])
Zpos_g <- c(8 + map$map_pos[map$V1 == 39])

indivlistgeno_A <- cbind.data.frame(indivlistgeno[,c(1:8)],indivlistgeno[,Apos_g])
indivlistgeno_Z <- cbind.data.frame(indivlistgeno[,c(1:8)],indivlistgeno[,Zpos_g])

save(indivlistgeno_A,file='working_files/intermediate_files/indivlistgeno_A.rdata')
save(indivlistgeno_Z,file='working_files/intermediate_files/indivlistgeno_Z.rdata')




