#script to model variance in allele frequency change over time - fixed sex ratio, autosomal version
#estimate sample values using real data
#Nancy Chen & Rose Driscoll
#Last updated: 5 May 2021

library(plyr)
library(tidyverse)

#get input files
load(file='../working_files/intermediate_files/indivlistgeno_A.rdata')
indivlistgeno <- indivlistgeno_A[,-c(8)]

#estimate values that vary with SNP
#n = number of genotyped individuals, x = sample allele frequency
#calculate allele freq
sampleFreq<-data.frame(Year=rep(c(1999:2013),each=33),Category=rep(c('nMt','xMt','nFt',
  'xFt','nMs','xMs','nFs','xFs','nMi','xMi','nFi','xFi','nMb','xMb','nFb','xFb',
  'xMt1-xMt','xFt1-xFt','xMs-xMt','xFs-xFt','xMi-xMt','xFi-xFt','xMb-xMt','xFb-xFt','xMdad',
  'xFdad','xMmom','xMfam','xFfam','xMmend','xFmend','xMfam-xMt','xFfam-xFt'),15),
  stringsAsFactors=FALSE)
# 'xMdad', 'xFdad', 'xMmom', and 'xFmom' are the frequencies with which males transmit an 
# allele to their sons, males transmit an allele to their daughters, females transmit an 
# allele to their sons, and females transmit an allele to their daughters respectively
# xFmom gets omitted for the Z as females never transmit their Z to their daughters, but is included for auto
sampleFreq<-rbind(data.frame(Year=rep(1998,4),Category=c('nMt','nFt','xMt','xFt'),
	stringsAsFactors=FALSE),sampleFreq)
SNPyr<-sampleFreq$Year
SNPcat<-sampleFreq$Category
igYear<-indivlistgeno$Year

#snp="SNP1"
for(snp in names(indivlistgeno)[8:length(indivlistgeno)])
{
  # number of genotyped chromosomes
	sampleFreq[SNPyr==1998 & SNPcat=='nMt',snp]<-
		2*sum(!is.na(indivlistgeno[igYear==1998&indivlistgeno$Sex==1,snp]))
	sampleFreq[SNPyr==1998 & SNPcat=='nFt',snp]<-
	  2*sum(!is.na(indivlistgeno[igYear==1998&indivlistgeno$Sex==2,snp]))
	
	#sampleFreq[SNPyr==1998 & SNPcat=='nt',snp]<-
	#  sampleFreq[SNPyr==1998 & SNPcat=='nMt',snp]+sampleFreq[SNPyr==1998 & SNPcat=='nFt',snp]

	#sampleFreq[SNPyr==1998 & SNPcat=='xt',snp]<-sum(indivlistgeno[igYear==1998,snp],
	#  na.rm=TRUE)/sampleFreq[SNPyr==1998 & SNPcat=='nt',snp]
	
	sampleFreq[SNPyr==1998 & SNPcat=='xMt',snp]<-sum(indivlistgeno[igYear==1998&indivlistgeno$Sex==1,snp],
		na.rm=TRUE)/sampleFreq[SNPyr==1998 & SNPcat=='nMt',snp]
	sampleFreq[SNPyr==1998 & SNPcat=='xFt',snp]<-sum(indivlistgeno[igYear==1998&indivlistgeno$Sex==2,snp],
	  na.rm=TRUE)/sampleFreq[SNPyr==1998 & SNPcat=='nFt',snp]

	for(year in c(1999:2013))
	{
		genoYr<-indivlistgeno[igYear==year,]
		
		sampleFreq[SNPyr==year & SNPcat=='nMt',snp]<-2*sum(!is.na(genoYr[genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFt',snp]<-2*sum(!is.na(genoYr[genoYr$Sex==2,snp]))
		
		#sampleFreq[SNPyr==year & SNPcat=='nt',snp]<-
		#  sampleFreq[SNPyr==year & SNPcat=='nMt',snp]+sampleFreq[SNPyr==year & SNPcat=='nFt',snp]

		#sampleFreq[SNPyr==year & SNPcat=='xt',snp]<-sum(genoYr[,snp],na.rm=TRUE)/
		#  sampleFreq[SNPyr==year & SNPcat=='nt',snp]
		
		sampleFreq[SNPyr==year & SNPcat=='xMt',snp]<-sum(genoYr[genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMt',snp]
		sampleFreq[SNPyr==year & SNPcat=='xFt',snp]<-sum(genoYr[genoYr$Sex==2,snp],na.rm=TRUE)/
		  sampleFreq[SNPyr==year & SNPcat=='nFt',snp]

		sampleFreq[SNPyr==year & SNPcat=='nMs',snp]<-
			2*sum(!is.na(genoYr[genoYr$Category=='survivor'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFs',snp]<-
		  2*sum(!is.na(genoYr[genoYr$Category=='survivor'&genoYr$Sex==2,snp]))
		
		sampleFreq[SNPyr==year & SNPcat=='xMs',snp]<-
			sum(genoYr[genoYr$Category=='survivor'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMs',snp]
		sampleFreq[SNPyr==year & SNPcat=='xFs',snp]<-
		  sum(genoYr[genoYr$Category=='survivor'&genoYr$Sex==2,snp],na.rm=TRUE)/
		  sampleFreq[SNPyr==year & SNPcat=='nFs',snp]
		
		sampleFreq[SNPyr==year & SNPcat=='nMi',snp]<-
			2*sum(!is.na(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFi',snp]<-
		  2*sum(!is.na(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==2,snp]))
		
		sampleFreq[SNPyr==year & SNPcat=='xMi',snp]<-
			ifelse(sampleFreq[SNPyr==year & SNPcat=='nMi',snp]==0,0,
			sum(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMi',snp])
		sampleFreq[SNPyr==year & SNPcat=='xFi',snp]<-
		  ifelse(sampleFreq[SNPyr==year & SNPcat=='nFi',snp]==0,0,
		  sum(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==2,snp],na.rm=TRUE)/
		  sampleFreq[SNPyr==year & SNPcat=='nFi',snp])
		
		sampleFreq[SNPyr==year & SNPcat=='nMb',snp]<-
			2*sum(!is.na(genoYr[genoYr$Category=='nestling'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFb',snp]<-
		  2*sum(!is.na(genoYr[genoYr$Category=='nestling'&genoYr$Sex==2,snp]))
		
		sampleFreq[SNPyr==year & SNPcat=='xMb',snp]<-
			sum(genoYr[genoYr$Category=='nestling'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMb',snp]
		sampleFreq[SNPyr==year & SNPcat=='xFb',snp]<-
		  sum(genoYr[genoYr$Category=='nestling'&genoYr$Sex==2,snp],na.rm=TRUE)/
		  sampleFreq[SNPyr==year & SNPcat=='nFb',snp]
	}
}

#calculate allele freq differences between each category and the year before
for(year in c(1999:2013))
{
	#allele freq differences
  #sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(4:252)]<-
  #  sampleFreq[SNPyr==year & SNPcat=='xt',c(4:252)]-
  #  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(4:252)]
  
	sampleFreq[SNPyr==year & SNPcat=='xMt1-xMt',c(3:length(sampleFreq))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMt',c(3:length(sampleFreq))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xMt',c(3:length(sampleFreq))]
	sampleFreq[SNPyr==year & SNPcat=='xFt1-xFt',c(3:length(sampleFreq))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFt',c(3:length(sampleFreq))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xFt',c(3:length(sampleFreq))]
	
	sampleFreq[SNPyr==year & SNPcat=='xMs-xMt',c(3:length(sampleFreq))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMs',c(3:length(sampleFreq))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xMt',c(3:length(sampleFreq))]
	sampleFreq[SNPyr==year & SNPcat=='xFs-xFt',c(3:length(sampleFreq))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFs',c(3:length(sampleFreq))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xFt',c(3:length(sampleFreq))]
	
	sampleFreq[SNPyr==year & SNPcat=='xMi-xMt',c(3:length(sampleFreq))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMi',c(3:length(sampleFreq))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xMt',c(3:length(sampleFreq))]
	sampleFreq[SNPyr==year & SNPcat=='xFi-xFt',c(3:length(sampleFreq))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFi',c(3:length(sampleFreq))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xFt',c(3:length(sampleFreq))]
	
	sampleFreq[SNPyr==year & SNPcat=='xMb-xMt',c(3:length(sampleFreq))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMb',c(3:length(sampleFreq))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xMt',c(3:length(sampleFreq))]
	sampleFreq[SNPyr==year & SNPcat=='xFb-xFt',c(3:length(sampleFreq))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFb',c(3:length(sampleFreq))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xFt',c(3:length(sampleFreq))]
}

#calculate variances and covariances
sampleVar<-data.frame(Year=rep(c(1999:2013),each=29),Category=rep(c('xMt1-xMt','xFt1-xFt',
  'xMs-xMt','xFs-xFt','xMi-xMt','xFi-xFt','xMb-xMt','xFb-xFt',
  'xMsxMi','xFsxFi','xMsxMb','xFsxFb','xMixMb','xFixFb',
  'xMsxFs','xMixFi','xMbxFb',
  'xMsxFi','xMsxFb','xMixFb',
  'xFsxMi','xFsxMb','xFixMb',
  'xMfam','xFfam','xMmend','xFmend','xMfam-xMt','xFfam-xFt'),15),stringsAsFactors=FALSE)
	
varYr<-sampleVar$Year
varCat<-sampleVar$Category

for(year in c(1999:2013))
{
	#var
  #sampleVar[varYr==year & varCat=='xt1-xt','avg']<-
  #  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(4:252)])^2)
  
	sampleVar[varYr==year & varCat=='xMt1-xMt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMt1-xMt',c(3:length(sampleFreq))])^2)
	sampleVar[varYr==year & varCat=='xFt1-xFt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFt1-xFt',c(3:length(sampleFreq))])^2)
		
	sampleVar[varYr==year & varCat=='xMs-xMt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xMt',c(3:length(sampleFreq))])^2)
	sampleVar[varYr==year & varCat=='xFs-xFt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xFt',c(3:length(sampleFreq))])^2)
	
	sampleVar[varYr==year & varCat=='xMi-xMt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xMt',c(3:length(sampleFreq))])^2)
	sampleVar[varYr==year & varCat=='xFi-xFt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xFt',c(3:length(sampleFreq))])^2)
	
	sampleVar[varYr==year & varCat=='xMb-xMt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xMt',c(3:length(sampleFreq))])^2)
	sampleVar[varYr==year & varCat=='xFb-xFt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xFt',c(3:length(sampleFreq))])^2)

	#covar
	sampleVar[varYr==year & varCat=='xMsxMi','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xMt',c(3:length(sampleFreq))])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xMt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xFsxFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xFt',c(3:length(sampleFreq))])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xFt',c(3:length(sampleFreq))]))
		
	sampleVar[varYr==year & varCat=='xMsxMb','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xMt',c(3:length(sampleFreq))])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xMt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xFsxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xFt',c(3:length(sampleFreq))])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xFt',c(3:length(sampleFreq))]))
		
	sampleVar[varYr==year & varCat=='xMixMb','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xMt',c(3:length(sampleFreq))])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xMt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xFixFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xFt',c(3:length(sampleFreq))])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xFt',c(3:length(sampleFreq))]))
	
	sampleVar[varYr==year & varCat=='xMsxFs','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xMt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xFt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xMixFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xMt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xFt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xMbxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xMt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xFt',c(3:length(sampleFreq))]))
	
	sampleVar[varYr==year & varCat=='xMsxFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xMt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xFt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xMsxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xMt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xFt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xMixFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xMt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xFt',c(3:length(sampleFreq))]))
	
	
	sampleVar[varYr==year & varCat=='xFsxMi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xFt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xMt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xFsxMb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xFt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xMt',c(3:length(sampleFreq))]))
	sampleVar[varYr==year & varCat=='xFixMb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xFt',c(3:length(sampleFreq))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xMt',c(3:length(sampleFreq))]))
	
}

#now calculate Mendelian noise
#get unique individuals
genoUnique<-indivlistgeno[!duplicated(indivlistgeno$Indiv),]

#get sample allele frequencies of parents
for(year in c(1999:2013))
{
	dadsofmales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==1,'Dad']
	dadsofmales<-data.frame(Indiv=dadsofmales[!is.na(dadsofmales)],stringsAsFactors=FALSE)
	dadsofmalesgeno<-merge(dadsofmales,genoUnique[,c(2,8:length(genoUnique))],by='Indiv',all.x=TRUE)
	
	momsofmales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==1,'Mom']
	momsofmales<-data.frame(Indiv=momsofmales[!is.na(momsofmales)],stringsAsFactors=FALSE)
	momsofmalesgeno<-merge(momsofmales,genoUnique[,c(2,8:length(genoUnique))],by='Indiv',all.x=TRUE)
	
	dadsoffemales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==2,'Dad']
	dadsoffemales<-data.frame(Indiv=dadsoffemales[!is.na(dadsoffemales)],stringsAsFactors=FALSE)
	dadsoffemalesgeno<-merge(dadsoffemales,genoUnique[,c(2,8:length(genoUnique))],by='Indiv',all.x=TRUE)
	
	#Unlike on the Z, moms *do* contribute to daughters on autosomes we must include momsoffemales
	momsoffemales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==2,'Mom']
	momsoffemales<-data.frame(Indiv=momsoffemales[!is.na(momsoffemales)],stringsAsFactors=FALSE)
	momsoffemalesgeno<-merge(momsoffemales,genoUnique[,c(2,8:length(genoUnique))],by='Indiv',all.x=TRUE)
	
	for(snp in names(indivlistgeno)[8:length(genoUnique)])
	{
		sampleFreq[SNPyr==year & SNPcat=='xMdad',snp]<-mean(dadsofmalesgeno[,snp],na.rm=TRUE)/2
		sampleFreq[SNPyr==year & SNPcat=='xMmom',snp]<-mean(momsofmalesgeno[,snp],na.rm=TRUE)/2
		
		sampleFreq[SNPyr==year & SNPcat=='xFdad',snp]<-mean(dadsoffemalesgeno[,snp],na.rm=TRUE)/2
		sampleFreq[SNPyr==year & SNPcat=='xFmom',snp]<-mean(momsoffemalesgeno[,snp],na.rm=TRUE)/2
	}
}

#calculate other terms
for(year in c(1999:2013))
{
	sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:length(sampleFreq))]<-0.5*
		(sampleFreq[SNPyr==year & SNPcat=='xMdad',c(3:length(sampleFreq))]+
		sampleFreq[SNPyr==year & SNPcat=='xMmom',c(3:length(sampleFreq))]) 
	  #here we multiply both mom and dad by 0.5 since each of their contributions makes up 1/2 of their son's genotype
	
	sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:length(sampleFreq))]<-0.5*
	  (sampleFreq[SNPyr==year & SNPcat=='xFdad',c(3:length(sampleFreq))]+
	  sampleFreq[SNPyr==year & SNPcat=='xFdad',c(3:length(sampleFreq))])
	  #here we multiply both mom and dad by 0.5 since each of their contributions makes up 1/2 of their daughter's genotype
		
	sampleFreq[SNPyr==year & SNPcat=='xMmend',c(3:length(sampleFreq))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMb',c(3:length(sampleFreq))]-
		sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:length(sampleFreq))]
	
	sampleFreq[SNPyr==year & SNPcat=='xFmend',c(3:length(sampleFreq))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFb',c(3:length(sampleFreq))]-
	  sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:length(sampleFreq))]
		
	sampleFreq[SNPyr==year & SNPcat=='xMfam-xMt',c(3:length(sampleFreq))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:length(sampleFreq))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xMt',c(3:length(sampleFreq))]
	
	sampleFreq[SNPyr==year & SNPcat=='xFfam-xFt',c(3:length(sampleFreq))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:length(sampleFreq))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xFt',c(3:length(sampleFreq))]
}

#calculate variances and covariances
for(year in c(1999:2013))
{
	sampleVar[varYr==year & varCat=='xMfam',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:length(sampleFreq))])^2)
	
	sampleVar[varYr==year & varCat=='xFfam',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:length(sampleFreq))])^2)
		
	sampleVar[varYr==year & varCat=='xMmend',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMmend',c(3:length(sampleFreq))])^2)
	
	sampleVar[varYr==year & varCat=='xFmend',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFmend',c(3:length(sampleFreq))])^2)
		
	sampleVar[varYr==year & varCat=='xMfam-xMt',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMfam-xMt',c(3:length(sampleFreq))])^2)
	
	sampleVar[varYr==year & varCat=='xFfam-xFt',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFfam-xFt',c(3:length(sampleFreq))])^2)
}

#get date script is run
today<-format(Sys.Date(),format="%d%b%Y")

#save output
save(sampleVar,file=paste("../working_files/intermediate_files/sampleVarA_FSR.rdata",sep=''))

