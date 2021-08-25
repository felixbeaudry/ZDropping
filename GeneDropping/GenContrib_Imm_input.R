#script to create ped files for estimating expected genetic contributions of immigrants
#Nancy Chen and Rose Driscoll
#Last updated: 3 April 2019

setwd('~/Documents/GitHub/ZDropping/ExpGenContrib')


library(plyr)

#consider contribution of cohorts (Yearly) or all immigrants (All)

#read in input files
pedigree<-read.table('working_files/pedigree.txt',header=FALSE,sep=' ',stringsAsFactors=FALSE)

#nestling data
indiv<-read.table('working_files/IndivDataUSFWS.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
cohortAll<-indiv[!is.na(indiv$NatalYear),1:2]

#assign immigrants from each year a different allele
#let's just use the year imm first appeared
imm<-pedigree[pedigree$Dad=='0',2:5]
imm$minObs<-laply(imm$Indiv,function(x) indiv[indiv$Indiv==x,'ImmCohort'])

#start in 1990
#all indiv assigned "0 0" genotypes
pedigree$allele1<-0
pedigree$allele2<-0

#assign "1 1" genotypes to immigrants from 1990 or before
pedigree[pedigree$Indiv %in% imm[,1],'allele1']<-1
pedigree[pedigree$Indiv %in% imm[,1],'allele2']<-1

#copy pedigree
pedigreeCopy<-pedigree

#immigrants in subsequent years (starting with 1990) given higher numbers. 
for (x in 1991:2013)
{
  # for males
	pedigree[pedigree$Indiv %in% imm[imm$minObs==x&imm$Sex==1,1],'allele1']<-2*(x-1990)
	pedigree[pedigree$Indiv %in% imm[imm$minObs==x&imm$Sex==1,1],'allele2']<-2*(x-1990)
	# for females
	pedigree[pedigree$Indiv %in% imm[imm$minObs==x&imm$Sex==2,1],'allele1']<-2*(x-1990)+1
	pedigree[pedigree$Indiv %in% imm[imm$minObs==x&imm$Sex==2,1],'allele2']<-2*(x-1990)+1
}
write.table(pedigree,file="working_files/intermediate_files/ImmContribYearly_malevsfemale.ped",quote=FALSE,sep=' ',row.names=FALSE,col.names=FALSE)

#immigrants in subsequent years (starting with 1991) given '2 2' genotypes if male and '3 3' genotypes if female
for (x in 1991:2013)
{
  # for males
	pedigreeCopy[pedigreeCopy$Indiv %in% imm[imm$minObs==x&imm$Sex==1,1],'allele1']<-2
	pedigreeCopy[pedigreeCopy$Indiv %in% imm[imm$minObs==x&imm$Sex==1,1],'allele2']<-2
	# for females
	pedigreeCopy[pedigreeCopy$Indiv %in% imm[imm$minObs==x&imm$Sex==2,1],'allele1']<-3
	pedigreeCopy[pedigreeCopy$Indiv %in% imm[imm$minObs==x&imm$Sex==2,1],'allele2']<-3
}
write.table(pedigreeCopy,file="working_files/intermediate_files/ImmContribAll_malevsfemale.ped",quote=FALSE,sep=' ',row.names=FALSE,col.names=FALSE)

#print out file with nestling cohorts
write.table(cohortAll,file="working_files/intermediate_files/allABSnestlings.txt",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

