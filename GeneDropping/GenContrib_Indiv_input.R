#script to create ped files for estimating expected genetic contributions of 926 breeders
#Nancy Chen and Rose Driscoll
#Last updated: 8 Sept 2021

library(plyr)
library(dplyr)

#read in input files
#1990-2013 breeders data
breeders <- read.table("working_files/926_breeders.txt", header = TRUE, sep = ' ', stringsAsFactors = FALSE)
#pedigree + genotype data
pedigree<-read.table("working_files/pedigree.txt",header=TRUE,sep=' ',stringsAsFactors=FALSE)

#generate ped files for each individual
for (ind in breeders$Indiv) {
	thisped<-pedigree
	
	#all indiv assigned "1 1" genotypes
	thisped$allele1<-1
	thisped$allele2<-1
		
	#except for focal indiv - given "2 2" genotypes
	thisped[thisped$Indiv == ind,'allele1']<-2
	thisped[thisped$Indiv == ind,'allele2']<-2
	
	#we also want to turn the focal indiv into a founder by changing its parents to 0
	thisped[thisped$Indiv == ind,'Dad']<-'0'
	thisped[thisped$Indiv == ind,'Mom']<-'0'
	
	write.table(thisped,file=paste("working_files/intermediate_files/IndivContrib_",ind,".ped",sep=''),quote=FALSE,sep=' ',row.names=FALSE,col.names=FALSE)
}

#get nestling data
indiv<-read.table('working_files/IndivDataUSFWS.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
cohortAll<-indiv[!is.na(indiv$NatalYear),1:2]

#print out file with nestling cohorts
write.table(cohortAll,file="working_files/intermediate_files/allABSnestlings.txt",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

