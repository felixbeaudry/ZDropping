#script to make input files for alleleFreqModel from plink and pedigree
#Felix Beaudry & Rose Driscoll & Nancy Chen
#Last updated: May 31 2022

library(tidyverse)
library(data.table)

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

###generate files for LD-pruned analyses
# LD pruning commands in PLINK
# plink --noweb --allow-no-sex --dog --nonfounders --file geno.anon --extract prune.in --recode --out geno.anon.pruned

# get pruned input files
ped_geno_LD <- read.table('genotypeFiles/geno.anon.pruned.ped', header = FALSE, sep = ' ', stringsAsFactors = FALSE)
map_LD <- read.table('genotypeFiles/geno.anon.pruned.map')

write.table(map_LD, 
            file = "alleleFreqModel/working_files/geno.pruned.map", append = F, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, qmethod = c("escape", "double"))

map_LD$map_pos <- seq(1, length(map_LD$V1))

# convert genotype data to one column per SNP
tmpped <- ped_geno_LD[ped_geno_LD$V7 != 0, ]
ped_geno_sumd_LD <- tmpped[, 2:5] 
# cycle through SNPs
for (x in seq(7, (6 + length(map_LD$V2) * 2), by = 2)) {
	thisped <- tmpped[, c(x, x+1)]
	thisped$geno <- rowSums(thisped[, 1:2]) - 2
	thisped[thisped$geno < 0, 'geno'] <- NA
	ped_geno_sumd_LD <- cbind(ped_geno_sumd_LD, thisped$geno)
}
names(ped_geno_sumd_LD)<-c('Indiv', 'Dad', 'Mom', 'Sex', paste("SNP",c(1:(length(map_LD$V2))), sep = ""))

# separate into autosomal and Z SNPs
ped_geno_sumd_A_LD <- ped_geno_sumd_LD[, c(1:4, 4 + map_LD$map_pos[map_LD$V1 %in% c(0:32)])]
ped_geno_sumd_Z_LD <- ped_geno_sumd_LD[, c(1:4, 4 + map_LD$map_pos[map_LD$V1 == 39])]

# combine
indivlistgenoA_LD <- left_join(indivlist, ped_geno_sumd_A_LD[,c(1, 5:length(ped_geno_sumd_A_LD))], by = c("Indiv" = "Indiv"))
indivlistgenoZ_LD <- left_join(indivlist, ped_geno_sumd_Z_LD[,c(1, 5:length(ped_geno_sumd_Z_LD))], by = c("Indiv" = "Indiv"))

# sex encoding: 1 = male, 2 = female
maleCode <- 1
femaleCode <- 2

# correct female genotypes on Z to indicate that they are heterozygous (1 instead of 2)
indivlistgenoZ_LD[indivlistgenoZ_LD$simSex == femaleCode, c(8:length(indivlistgenoZ_LD))] <- indivlistgenoZ_LD[indivlistgenoZ_LD$simSex == femaleCode, c(8:length(indivlistgenoZ_LD))] / 2

# export files
save(indivlistgenoA_LD, file = 'alleleFreqModel/working_files/intermediate_files/indivlistgeno_Apruned.rdata')
save(indivlistgenoZ_LD, file = 'alleleFreqModel/working_files/intermediate_files/indivlistgeno_Zpruned.rdata')