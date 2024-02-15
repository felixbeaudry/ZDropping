#script to make input files for alleleFreqModel from plink and pedigree
#Felix Beaudry & Rose Driscoll & Nancy Chen
#Last updated: May 31 2022

library(tidyverse)
library(data.table)
'%ni%' <- Negate('%in%')

####top####

ped_geno <- read.table('genotypeFiles/geno.anon.ped', header = FALSE, sep = ' ', stringsAsFactors = FALSE)
ped <- ped_geno[,c(1:6)]

write.table(ped, 
            file = "alleleFreqModel/working_files/pedigree.txt", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, qmethod = c("escape", "double"))

write.table(ped, 
            file = "GeneDropping/working_files/pedigree.txt", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, qmethod = c("escape", "double"))

map<-read.table('genotypeFiles/geno.anon.map')
write.table(map, 
            file = "alleleFreqModel/working_files/geno.map", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, qmethod = c("escape", "double"))
map$map_pos <- seq(1, length(map$V1))

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

####make pedgenos for allele freq models####
indivlist <- fread('genotypeFiles/indivlist.txt')

# We use simulated sexes for unsexed individuals in the pedigree in our allele frequency models.
# For reproducibility, the simulated sex data that we use is preserved in column `simSex` in indivlist.txt
# Below, we show but do not run the code that was used to simulate sexes.

#check for unsexed indivs & assign them a sex
#unsexed_indivs <- ped$Indiv[ped$Sex==0]
#simulated_sexes <- sample(x = c(1,2), size = length(unsexed_indivs), prob = c(0.5,0.5), replace = TRUE)
#ped$simSex <- ped$Sex
#ped$simSex[ped$simSex==0] <- simulated_sexes

#indivlist <- left_join(indivlist,ped,by=c("USFWS"="Indiv"))
#names(indivlist)<-c('Year','Indiv','Category','Genotyped','Mom','Dad','Sex','simSex')

# convert genotype data to one column per SNP
tmpped <- ped_geno[ped_geno$V7 != 0, ]
ped_geno_sumd <- tmpped[, 2:5] 
# cycle through SNPs
for (x in seq(7, (6 + length(map$V2) * 2), by = 2)) {
	thisped <- tmpped[, c(x, x+1)]
	thisped$geno <- rowSums(thisped[, 1:2]) - 2
	thisped[thisped$geno < 0, 'geno'] <- NA
	ped_geno_sumd <- cbind(ped_geno_sumd, thisped$geno)
}
names(ped_geno_sumd)<-c('Indiv', 'Dad', 'Mom', 'Sex', paste("SNP",c(1:(length(map$V2))), sep = ""))

# separate into autosomal and Z SNPs
ped_geno_sumd_A <- ped_geno_sumd[, c(1:4, 4 + map$map_pos[map$V1 %in% c(0:32)])]
ped_geno_sumd_Z <- ped_geno_sumd[, c(1:4, 4 + map$map_pos[map$V1 == 39])]

# combine
indivlistgenoA <- left_join(indivlist, ped_geno_sumd_A[,c(1, 5:length(ped_geno_sumd_A))], by = c("Indiv" = "Indiv"))
indivlistgenoZ <- left_join(indivlist, ped_geno_sumd_Z[,c(1, 5:length(ped_geno_sumd_Z))], by = c("Indiv" = "Indiv"))

# sex encoding: 1 = male, 2 = female
maleCode <- 1
femaleCode <- 2

# correct female genotypes on Z to indicate that they are heterozygous (1 instead of 2)
indivlistgenoZ[indivlistgenoZ$simSex == femaleCode, c(8:length(indivlistgenoZ))] <- indivlistgenoZ[indivlistgenoZ$simSex == femaleCode, c(8:length(indivlistgenoZ))] / 2

save(indivlistgeno_A,file='alleleFreqModel/working_files/intermediate_files/indivlistgeno_A.rdata')
save(indivlistgeno_Z,file='alleleFreqModel/working_files/intermediate_files/indivlistgeno_Z.rdata')
