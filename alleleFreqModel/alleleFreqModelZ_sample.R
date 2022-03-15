#script to model variance in allele frequency change over time
#estimate sample values using real data for Z loci
#Nancy Chen, Rose Driscoll & Felix Beaudry
#tested on R v.4.1.2

library(tidyverse) #v.1.3.1
`%notin%` <- Negate(`%in%`)

####get & make starting data.frames####

#get input files
#setwd('~/Documents/Github/ZDropping/alleleFreqModel')
load(file='working_files/intermediate_files/indivlistgeno_Z.rdata')

ldprune.SNP <- read.table('ldprune.Z.SNP.list', header = FALSE, sep = "", dec = ".")
ldprune.cols <- c(names(indivlistgeno_Z)[c(1:7)],ldprune.SNP$V1)
indivlistgeno <- indivlistgeno_Z[,ldprune.cols]

snp_length <- length(indivlistgeno)-7

####allele frequencies####
#n = number of genotyped individuals, x = sample allele frequency
#calculate allele freq
sampleFreq<-data.frame(Year=rep(c(1999:2013),each=30),Category=rep(c(
  'nt','xt','xt1-xt',
  'nMs','xMs','nFs','xFs', 'xMs-xt','xFs-xt',
  'nMi','xMi','nFi','xFi', 'xMi-xt','xFi-xt',
  'nMb','xMb','nFb','xFb', 'xMb-xt','xFb-xt',
  'xMdad','xFdad','xMmom','xMfam','xFfam','xMmend','xFmend','xMfam-xt','xFfam-xt'),15),
  stringsAsFactors=FALSE)

# 'xMdad', 'xFdad', 'xMmom', and 'xFmom' are the frequencies with which males transmit an 
# allele to their sons, males transmit an allele to their daughters, females transmit an 
# allele to their sons, and females transmit an allele to their daughters respectively
# Omit xFmom for the Z as females never transmit their Z to their daughters

sampleFreq<-rbind.data.frame(data.frame(Year=rep(1998,2),Category=c('nt','xt'),
	stringsAsFactors=FALSE),sampleFreq)

for(snp in names(indivlistgeno)[8:(snp_length+7)]){
  sampleFreq <- sampleFreq %>% add_column("{snp}" := NA)
}


#make loop variables
SNPyr<-sampleFreq$Year
SNPcat<-sampleFreq$Category
igYear<-indivlistgeno$Year

#for each SNP, calculate number of genotyped chromosomes and use that to calculate allele frequency in each year
for(snp in names(indivlistgeno)[8:(snp_length+7)]){

  # number of genotyped chromosomes
	sampleFreq[SNPyr==1998 & SNPcat=='nt',snp]<-
		2*sum(!is.na(indivlistgeno[igYear==1998&indivlistgeno$Sex==1,snp])) +
	  sum(!is.na(indivlistgeno[igYear==1998&indivlistgeno$Sex==2,snp]))
	
	#allele frequency
	sampleFreq[SNPyr==1998 & SNPcat=='xt',snp]<-
	  sum(indivlistgeno[igYear==1998,snp],na.rm=TRUE)/
	  sampleFreq[SNPyr==1998 & SNPcat=='nt',snp]
	
	for(year in c(1999:2013)){
		genoYr<-indivlistgeno[igYear==year,]
		
		# number of genotyped chromosomes
		sampleFreq[SNPyr==year & SNPcat=='nt',snp]<-
		  2*sum(!is.na(genoYr[genoYr$Sex==1,snp])) + sum(!is.na(genoYr[genoYr$Sex==2,snp]))
		sampleFreq[SNPyr==year & SNPcat=='xt',snp]<-sum(genoYr[,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nt',snp]

		sampleFreq[SNPyr==year & SNPcat=='nMs',snp]<-
			2*sum(!is.na(genoYr[genoYr$Category=='survivor'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFs',snp]<-
		  sum(!is.na(genoYr[genoYr$Category=='survivor'&genoYr$Sex==2,snp]))
		
		sampleFreq[SNPyr==year & SNPcat=='nMi',snp]<-
		  2*sum(!is.na(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFi',snp]<-
		  sum(!is.na(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==2,snp]))
		
		sampleFreq[SNPyr==year & SNPcat=='nMb',snp]<-
		  2*sum(!is.na(genoYr[genoYr$Category=='nestling'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFb',snp]<-
		  sum(!is.na(genoYr[genoYr$Category=='nestling'&genoYr$Sex==2,snp]))
		
		#allele frequency
		sampleFreq[SNPyr==year & SNPcat=='xMs',snp]<-
			sum(genoYr[genoYr$Category=='survivor'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMs',snp]
		sampleFreq[SNPyr==year & SNPcat=='xFs',snp]<-
		  sum(genoYr[genoYr$Category=='survivor'&genoYr$Sex==2,snp],na.rm=TRUE)/
		  sampleFreq[SNPyr==year & SNPcat=='nFs',snp]
		
		sampleFreq[SNPyr==year & SNPcat=='xMi',snp]<-
			ifelse(sampleFreq[SNPyr==year & SNPcat=='nMi',snp]==0,0,
			sum(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMi',snp])
		sampleFreq[SNPyr==year & SNPcat=='xFi',snp]<-
		  ifelse(sampleFreq[SNPyr==year & SNPcat=='nFi',snp]==0,0,
		         sum(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==2,snp],na.rm=TRUE)/
		           sampleFreq[SNPyr==year & SNPcat=='nFi',snp])
		
		sampleFreq[SNPyr==year & SNPcat=='xMb',snp]<-
			sum(genoYr[genoYr$Category=='nestling'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMb',snp]
		sampleFreq[SNPyr==year & SNPcat=='xFb',snp]<-
		  sum(genoYr[genoYr$Category=='nestling'&genoYr$Sex==2,snp],na.rm=TRUE)/
		  sampleFreq[SNPyr==year & SNPcat=='nFb',snp]
	}
}

####change in allele frequencies####
#calculate allele freq differences between each category and the year before
for(year in c(1999:2013)){
	#allele freq differences
	sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(3:(snp_length+2))]<-
		sampleFreq[SNPyr==year & SNPcat=='xt',c(3:(snp_length+2))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
	
	sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:(snp_length+2))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMs',c(3:(snp_length+2))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
	sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:(snp_length+2))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFs',c(3:(snp_length+2))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
	
	sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:(snp_length+2))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMi',c(3:(snp_length+2))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
	sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:(snp_length+2))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFi',c(3:(snp_length+2))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
	
	sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:(snp_length+2))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMb',c(3:(snp_length+2))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
	sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:(snp_length+2))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFb',c(3:(snp_length+2))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
}

#calculate variances and covariances
sampleVar<-data.frame(Year=rep(c(1999:2013),each=28),Category=rep(c(
  'xt1-xt',
  'xMs-xt','xFs-xt',
  'xMi-xt','xFi-xt',
  'xMb-xt','xFb-xt',
  'xMsxMi','xFsxFi','xMsxMb','xFsxFb','xMixMb','xFixFb',
  'xMsxFs','xMixFi','xMbxFb',
  'xMsxFi','xMsxFb','xMixFb',
  'xFsxMi','xFsxMb','xFixMb',
  'xMfam','xFfam','xMmend','xFmend','xMfam-xt','xFfam-xt'),15),stringsAsFactors=FALSE)
	
varYr<-sampleVar$Year
varCat<-sampleVar$Category

for(year in c(1999:2013)){
	#var
	sampleVar[varYr==year & varCat=='xt1-xt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(3:(snp_length+2))])^2)
		
	sampleVar[varYr==year & varCat=='xMs-xt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:(snp_length+2))])^2)
	sampleVar[varYr==year & varCat=='xFs-xt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:(snp_length+2))])^2)
	
	sampleVar[varYr==year & varCat=='xMi-xt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:(snp_length+2))])^2)
	sampleVar[varYr==year & varCat=='xFi-xt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:(snp_length+2))])^2)
	
	sampleVar[varYr==year & varCat=='xMb-xt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:(snp_length+2))])^2)
	sampleVar[varYr==year & varCat=='xFb-xt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:(snp_length+2))])^2)

	#covar
	sampleVar[varYr==year & varCat=='xMsxMi','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:(snp_length+2))])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xFsxFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:(snp_length+2))])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:(snp_length+2))]))
		
	sampleVar[varYr==year & varCat=='xMsxMb','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:(snp_length+2))])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xFsxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:(snp_length+2))])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:(snp_length+2))]))
		
	sampleVar[varYr==year & varCat=='xMixMb','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:(snp_length+2))])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xFixFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:(snp_length+2))])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:(snp_length+2))]))
	
	sampleVar[varYr==year & varCat=='xMsxFs','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xMixFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xMbxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:(snp_length+2))]))
	
	sampleVar[varYr==year & varCat=='xMsxFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xMsxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xMixFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:(snp_length+2))]))
	
	sampleVar[varYr==year & varCat=='xFsxMi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xFsxMb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:(snp_length+2))]))
	sampleVar[varYr==year & varCat=='xFixMb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:(snp_length+2))])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:(snp_length+2))]))
	
}

####Mendelian noise####
#get unique individuals
genoUnique<-indivlistgeno[!duplicated(indivlistgeno$Indiv),]
#names(genoUnique)[2] <- "USFWS"

#year=1999
#get sample allele frequencies of parents
for(year in c(1999:2013)){
  
	dadsofmales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==1,'Dad']
	dadsofmalesgeno<-genoUnique[,c(2,8:(snp_length+7))] %>% filter(Indiv %in% dadsofmales)
	
	momsofmales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==1,'Mom']
	momsofmalesgeno<-genoUnique[,c(2,8:(snp_length+7))] %>% filter(Indiv %in% momsofmales)
	
	dadsoffemales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==2,'Dad']
	dadsoffemalesgeno<-genoUnique[,c(2,8:(snp_length+7))] %>% filter(Indiv %in% dadsoffemales)
	
	#Moms don't contribute to daughters on Z so skip these
	#momsoffemales<-indivlist[indivlist$Year==year & indivlist$Category=='nestling' & indivlist$Sex==2,'Mom']
	#momsoffemales<-data.frame(USFWS=momsoffemales[!is.na(momsoffemales)],stringsAsFactors=FALSE)
	#momsoffemalesgeno<-merge(momsoffemales,genoUnique[,c(1,8:256)],by='USFWS',all.x=TRUE)
	
	for(snp in names(indivlistgeno)[8:(snp_length+7)]){
		sampleFreq[SNPyr==year & SNPcat=='xMdad',snp]<-mean(unlist(dadsofmalesgeno[,snp]),na.rm=TRUE)/2
		sampleFreq[SNPyr==year & SNPcat=='xMmom',snp]<-mean(unlist(momsofmalesgeno[,snp]),na.rm=TRUE)
		
		sampleFreq[SNPyr==year & SNPcat=='xFdad',snp]<-mean(unlist(dadsoffemalesgeno[,snp]),na.rm=TRUE)/2
		#sampleFreq[SNPyr==year & SNPcat=='xFmom',snp]<-mean(momsoffemalesgeno[,snp],na.rm=TRUE)/2
	}
}

#calculate other terms
for(year in c(1999:2013)){
	sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:(snp_length+2))]<-0.5*
		(sampleFreq[SNPyr==year & SNPcat=='xMdad',c(3:(snp_length+2))]+
		sampleFreq[SNPyr==year & SNPcat=='xMmom',c(3:(snp_length+2))]) 

	sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:(snp_length+2))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFdad',c(3:(snp_length+2))]

	sampleFreq[SNPyr==year & SNPcat=='xMmend',c(3:(snp_length+2))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMb',c(3:(snp_length+2))]-
		sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:(snp_length+2))]
	
	sampleFreq[SNPyr==year & SNPcat=='xFmend',c(3:(snp_length+2))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFb',c(3:(snp_length+2))]-
	  sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:(snp_length+2))]
		
	sampleFreq[SNPyr==year & SNPcat=='xMfam-xt',c(3:(snp_length+2))]<-
		sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:(snp_length+2))]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
	
	sampleFreq[SNPyr==year & SNPcat=='xFfam-xt',c(3:(snp_length+2))]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:(snp_length+2))]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:(snp_length+2))]
}

#calculate variances and covariances
for(year in c(1999:2013)){
	sampleVar[varYr==year & varCat=='xMfam',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:(snp_length+2))])^2)
	
	sampleVar[varYr==year & varCat=='xFfam',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:(snp_length+2))])^2)
		
	sampleVar[varYr==year & varCat=='xMmend',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMmend',c(3:(snp_length+2))])^2)
	
	sampleVar[varYr==year & varCat=='xFmend',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFmend',c(3:(snp_length+2))])^2)
		
	sampleVar[varYr==year & varCat=='xMfam-xt',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMfam-xt',c(3:(snp_length+2))])^2)
	
	sampleVar[varYr==year & varCat=='xFfam-xt',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFfam-xt',c(3:(snp_length+2))])^2)
}


#save output
save(sampleVar,file=paste("working_files/intermediate_files/sampleVar_Z.ldprune.rdata",sep=''))


