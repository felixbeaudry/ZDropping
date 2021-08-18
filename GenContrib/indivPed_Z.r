#script to create ped files for estimating expected genetic contributions of founder individuals
#Nancy Chen and Rose Driscoll
#Last updated: 1 Apr 2019

library(plyr)

#read in input files
#pedigree + genotype data
ped<-read.table('FSJpedgeno_Zsexlinked.ped',header=TRUE,sep=' ',stringsAsFactors=FALSE)
pedigree<-ped[,1:6]
names(pedigree)<-c('Fam','Indiv','Dad','Mom','Sex','Pheno')

#get founders
founders<-data.frame(USFWS=pedigree[pedigree$Dad=='0','Indiv'],stringsAsFactors=FALSE)

#only need to run sims for founders with kids
founders$kids<-founders$USFWS %in% c(pedigree$Dad,pedigree$Mom)	

#generate ped files for each individual
for (ind in founders[founders$kids,'USFWS']) {
	thisped<-pedigree
	
	#all indiv assigned "1 1" genotypes
	thisped$allele1<-1
	thisped$allele2<-1
		
	#except for focal indiv - given "2 2" genotypes
	thisped[thisped$Indiv == ind,'allele1']<-2
	thisped[thisped$Indiv == ind,'allele2']<-2
	
	write.table(thisped,file=paste("IndivContrib_",ind,".ped",sep=''),quote=FALSE,sep=' ',row.names=FALSE,col.names=FALSE)
}

#get nestling data
indiv<-read.table('IndivDataUSFWS.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
cohortAll<-indiv[!is.na(indiv$NatalYear),1:2]

#print out file with nestling cohorts
write.table(cohortAll,file="allABSnestlings.txt",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

