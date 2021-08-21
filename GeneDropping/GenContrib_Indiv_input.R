#script to create ped files for estimating expected genetic contributions of founder individuals
#Nancy Chen and Rose Driscoll
#Last updated: 1 Apr 2019

library(plyr)

#read in input files
pedigree<-read.table('working_files/pedigree.txt',header=FALSE,sep=' ',stringsAsFactors=FALSE)


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
	
	write.table(thisped,file=paste("working_files/intermediate_files/IndivContrib_",ind,".ped",sep=''),quote=FALSE,sep=' ',row.names=FALSE,col.names=FALSE)
}

#get nestling data
indiv<-read.table('working_files/IndivDataUSFWS.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
cohortAll<-indiv[!is.na(indiv$NatalYear),1:2]

#print out file with nestling cohorts
write.table(cohortAll,file="working_files/intermediate_files/allABSnestlings.txt",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

