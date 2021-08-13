#script to model variance in allele frequency change over time
#estimate sample values using real data for Z loci
#Nancy Chen, Rose Driscoll & Felix Beaudry
#Last updated: 7 July 2021

library(plyr)


####get & make starting data.frames####
#list of indiv in each category each year
load("simindivFIXmin2obs.rdata")

#genotype data
ped<-read.table('FSJpedgeno_Zsexlinked.ped',header=FALSE,sep=' ',stringsAsFactors=FALSE)

#subset pedigree for genotyped individuals
Genotyped_official <- unique(simindivFIXmin2obs$USFWS[simindivFIXmin2obs$genotyped == "Y"])
tmpped<-ped[ped$V2 %in% Genotyped_official,]
pedgeno<-tmpped[,2:5]

#convert genotype data to one column per SNP
for (x in seq(7,504,by=2)) { #cycle through SNPs
	thisped<-tmpped[,c(x,x+1)]
	thisped$geno<-rowSums(thisped[,1:2])-2
	thisped[thisped$geno<0,'geno']<-NA
	pedgeno<-cbind(pedgeno,thisped$geno)
}
names(pedgeno)<-c('Indiv','Dad','Mom','Sex',c(1:249))

#add genotype info to pedigree
indivlistgeno<-merge(simindivFIXmin2obs[,1:4],pedgeno[,c(1:4,5:253)],by.x='USFWS',by.y='Indiv')


#correct female genotypes to indicate that they are heterozygous (1 instead of 2)
indivlistgeno[indivlistgeno[,7]==2,8:256] <- indivlistgeno[indivlistgeno[,7]==2,8:256]/2

#save genotype file for simulations
save(indivlistgeno,file='indivlistgenoZ.rdata')

#add sex to indivlist
indivlist <- merge(simindivFIXmin2obs,pedgeno[c(1,4)],by.x="USFWS",by.y="Indiv")
names(indivlist)<-c('Indiv','Year','Category','Genotyped','Mom','Dad','Sex')

####values constant across SNPs####
samplePars<-data.frame(Year=c(1998:2013),stringsAsFactors=FALSE)

samplePars[samplePars$Year==1998,'Nt']<- 
  2*length(indivlist[indivlist$Year==1998&indivlist$Category=='survivor'&indivlist$Sex==1,1]) + 
  2*length(indivlist[indivlist$Year==1998&indivlist$Category=='nestling'&indivlist$Sex==1,1]) +
  2*length(indivlist[indivlist$Year==1998&indivlist$Category=='immigrant'&indivlist$Sex==1,1]) +
  length(indivlist[indivlist$Year==1998&indivlist$Category=='survivor'&indivlist$Sex==2,1]) + 
  length(indivlist[indivlist$Year==1998&indivlist$Category=='immigrant'&indivlist$Sex==2,1]) + 
  length(indivlist[indivlist$Year==1998&indivlist$Category=='nestling'&indivlist$Sex==2,1]) 

#calculate chromosomes by demographic category for each year
for(yr in c(1999:2013)){
  
	samYr<-samplePars$Year
	indivYr<-indivlist[which(indivlist$Year==yr),]
	
	#number of inds each yr in each category (total, survivors, immigrants, nestlings of each sex)
	samplePars[samYr==yr,'NMs']<-2*length(indivYr[indivYr$Category=='survivor'&indivYr$Sex==1,1])
	samplePars[samYr==yr,'NFs']<-length(indivYr[indivYr$Category=='survivor'&indivYr$Sex==2,1])
	samplePars[samYr==yr,'NMi']<-2*length(indivYr[indivYr$Category=='immigrant'&indivYr$Sex==1,1])
	samplePars[samYr==yr,'NFi']<-length(indivYr[indivYr$Category=='immigrant'&indivYr$Sex==2,1])
	samplePars[samYr==yr,'NMb']<-2*length(indivYr[indivYr$Category=='nestling'&indivYr$Sex==1,1])
	samplePars[samYr==yr,'NFb']<-length(indivYr[indivYr$Category=='nestling'&indivYr$Sex==2,1])
	
	samplePars[samYr==yr,'Nt']<- samplePars[samYr==yr,'NFb'] + samplePars[samYr==yr,'NMb'] + samplePars[samYr==yr,'NFi'] + 	samplePars[samYr==yr,'NMi'] + samplePars[samYr==yr,'NMs'] +  samplePars[samYr==yr,'NFs'] 
	
	#proportion of chromosomes each yr in each category
	samplePars[samYr==yr,'propMS']<-samplePars[samYr==yr,'NMs']/samplePars[samYr==yr,'Nt']
	samplePars[samYr==yr,'propFS']<-samplePars[samYr==yr,'NFs']/samplePars[samYr==yr,'Nt']
	samplePars[samYr==yr,'propMI']<-samplePars[samYr==yr,'NMi']/samplePars[samYr==yr,'Nt']
	samplePars[samYr==yr,'propFI']<-samplePars[samYr==yr,'NFi']/samplePars[samYr==yr,'Nt']
	samplePars[samYr==yr,'propMB']<-samplePars[samYr==yr,'NMb']/samplePars[samYr==yr,'Nt']
	samplePars[samYr==yr,'propFB']<-samplePars[samYr==yr,'NFb']/samplePars[samYr==yr,'Nt']
}

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

sampleFreq<-rbind(data.frame(Year=rep(1998,2),Category=c('nt','xt'),
	stringsAsFactors=FALSE),sampleFreq)

#make loop variables
SNPyr<-sampleFreq$Year
SNPcat<-sampleFreq$Category
igYear<-indivlistgeno$Year

#for each SNP, calculate number of genotyped chromosomes and use that to calculate allele frequency in each year
for(snp in names(indivlistgeno)[8:256]){

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
			2*sum(!is.na(genoYr[genoYr$category=='survivor'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFs',snp]<-
		  sum(!is.na(genoYr[genoYr$category=='survivor'&genoYr$Sex==2,snp]))
		
		sampleFreq[SNPyr==year & SNPcat=='nMi',snp]<-
		  2*sum(!is.na(genoYr[genoYr$category=='immigrant'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFi',snp]<-
		  sum(!is.na(genoYr[genoYr$category=='immigrant'&genoYr$Sex==2,snp]))
		
		sampleFreq[SNPyr==year & SNPcat=='nMb',snp]<-
		  2*sum(!is.na(genoYr[genoYr$category=='nestling'&genoYr$Sex==1,snp]))
		sampleFreq[SNPyr==year & SNPcat=='nFb',snp]<-
		  sum(!is.na(genoYr[genoYr$category=='nestling'&genoYr$Sex==2,snp]))
		
		#allele frequency
		sampleFreq[SNPyr==year & SNPcat=='xMs',snp]<-
			sum(genoYr[genoYr$category=='survivor'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMs',snp]
		sampleFreq[SNPyr==year & SNPcat=='xFs',snp]<-
		  sum(genoYr[genoYr$category=='survivor'&genoYr$Sex==2,snp],na.rm=TRUE)/
		  sampleFreq[SNPyr==year & SNPcat=='nFs',snp]
		
		sampleFreq[SNPyr==year & SNPcat=='xMi',snp]<-
			ifelse(sampleFreq[SNPyr==year & SNPcat=='nMi',snp]==0,0,
			sum(genoYr[genoYr$category=='immigrant'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMi',snp])
		sampleFreq[SNPyr==year & SNPcat=='xFi',snp]<-
		  ifelse(sampleFreq[SNPyr==year & SNPcat=='nFi',snp]==0,0,
		         sum(genoYr[genoYr$category=='immigrant'&genoYr$Sex==2,snp],na.rm=TRUE)/
		           sampleFreq[SNPyr==year & SNPcat=='nFi',snp])
		
		sampleFreq[SNPyr==year & SNPcat=='xMb',snp]<-
			sum(genoYr[genoYr$category=='nestling'&genoYr$Sex==1,snp],na.rm=TRUE)/
			sampleFreq[SNPyr==year & SNPcat=='nMb',snp]
		sampleFreq[SNPyr==year & SNPcat=='xFb',snp]<-
		  sum(genoYr[genoYr$category=='nestling'&genoYr$Sex==2,snp],na.rm=TRUE)/
		  sampleFreq[SNPyr==year & SNPcat=='nFb',snp]
	}
}

####change in allele frequencies####
#calculate allele freq differences between each category and the year before
for(year in c(1999:2013)){
	#allele freq differences
	sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(3:251)]<-
		sampleFreq[SNPyr==year & SNPcat=='xt',c(3:251)]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
	sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(3:251)]<-
	  sampleFreq[SNPyr==year & SNPcat=='xt',c(3:251)]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
	
	sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:251)]<-
		sampleFreq[SNPyr==year & SNPcat=='xMs',c(3:251)]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
	sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:251)]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFs',c(3:251)]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
	
	sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:251)]<-
		sampleFreq[SNPyr==year & SNPcat=='xMi',c(3:251)]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
	sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:251)]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFi',c(3:251)]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
	
	sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:251)]<-
		sampleFreq[SNPyr==year & SNPcat=='xMb',c(3:251)]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
	sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:251)]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFb',c(3:251)]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
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
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(3:251)])^2)
		
	sampleVar[varYr==year & varCat=='xMs-xt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:251)])^2)
	sampleVar[varYr==year & varCat=='xFs-xt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:251)])^2)
	
	sampleVar[varYr==year & varCat=='xMi-xt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:251)])^2)
	sampleVar[varYr==year & varCat=='xFi-xt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:251)])^2)
	
	sampleVar[varYr==year & varCat=='xMb-xt','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:251)])^2)
	sampleVar[varYr==year & varCat=='xFb-xt','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:251)])^2)

	#covar
	sampleVar[varYr==year & varCat=='xMsxMi','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:251)])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xFsxFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:251)])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:251)]))
		
	sampleVar[varYr==year & varCat=='xMsxMb','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:251)])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xFsxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:251)])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:251)]))
		
	sampleVar[varYr==year & varCat=='xMixMb','avg']<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:251)])*
		as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xFixFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:251)])*
	  as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:251)]))
	
	sampleVar[varYr==year & varCat=='xMsxFs','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xMixFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xMbxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:251)]))
	
	sampleVar[varYr==year & varCat=='xMsxFi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xMsxFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xMixFb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:251)]))
	
	sampleVar[varYr==year & varCat=='xFsxMi','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xFsxMb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:251)]))
	sampleVar[varYr==year & varCat=='xFixMb','avg']<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:251)])*
	         as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:251)]))
	
}

####Mendelian noise####
#get unique individuals
genoUnique<-indivlistgeno[!duplicated(indivlistgeno$USFWS),]

#get sample allele frequencies of parents
for(year in c(1999:2013)){
	dadsofmales<-indivlist[indivlist$Year==year & indivlist$Category=='nestling' & indivlist$Sex==1,'Dad']
	dadsofmales<-data.frame(USFWS=dadsofmales[!is.na(dadsofmales)],stringsAsFactors=FALSE)
	dadsofmalesgeno<-merge(dadsofmales,genoUnique[,c(1,8:256)],by='USFWS',all.x=TRUE)
	
	momsofmales<-indivlist[indivlist$Year==year & indivlist$Category=='nestling' & indivlist$Sex==1,'Mom']
	momsofmales<-data.frame(USFWS=momsofmales[!is.na(momsofmales)],stringsAsFactors=FALSE)
	momsofmalesgeno<-merge(momsofmales,genoUnique[,c(1,8:256)],by='USFWS',all.x=TRUE)
	
	dadsoffemales<-indivlist[indivlist$Year==year & indivlist$Category=='nestling' & indivlist$Sex==2,'Dad']
	dadsoffemales<-data.frame(USFWS=dadsoffemales[!is.na(dadsoffemales)],stringsAsFactors=FALSE)
	dadsoffemalesgeno<-merge(dadsoffemales,genoUnique[,c(1,8:256)],by='USFWS',all.x=TRUE)
	
	#Moms don't contribute to daughters on Z so skip these
	#momsoffemales<-indivlist[indivlist$Year==year & indivlist$Category=='nestling' & indivlist$Sex==2,'Mom']
	#momsoffemales<-data.frame(USFWS=momsoffemales[!is.na(momsoffemales)],stringsAsFactors=FALSE)
	#momsoffemalesgeno<-merge(momsoffemales,genoUnique[,c(1,8:256)],by='USFWS',all.x=TRUE)
	
	for(snp in names(indivlistgeno)[8:256])
	{
		sampleFreq[SNPyr==year & SNPcat=='xMdad',snp]<-mean(dadsofmalesgeno[,snp],na.rm=TRUE)/2
		sampleFreq[SNPyr==year & SNPcat=='xMmom',snp]<-mean(momsofmalesgeno[,snp],na.rm=TRUE)
		
		sampleFreq[SNPyr==year & SNPcat=='xFdad',snp]<-mean(dadsoffemalesgeno[,snp],na.rm=TRUE)/2
		#sampleFreq[SNPyr==year & SNPcat=='xFmom',snp]<-mean(momsoffemalesgeno[,snp],na.rm=TRUE)/2
	}
}

#calculate other terms
for(year in c(1999:2013)){
	sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:251)]<-0.5*
		(sampleFreq[SNPyr==year & SNPcat=='xMdad',c(3:251)]+
		sampleFreq[SNPyr==year & SNPcat=='xMmom',c(3:251)]) 
	  #here we multiply mom by 0.5 since her contribution makes up only 1/2 of her son's genotype
	
	sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:251)]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFdad',c(3:251)]
	  #here we don't multiply dad by 0.5 since his contribution makes up 100% of his daughter's genotype
		
	sampleFreq[SNPyr==year & SNPcat=='xMmend',c(3:251)]<-
		sampleFreq[SNPyr==year & SNPcat=='xMb',c(3:251)]-
		sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:251)]
	
	sampleFreq[SNPyr==year & SNPcat=='xFmend',c(3:251)]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFb',c(3:251)]-
	  sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:251)]
		
	sampleFreq[SNPyr==year & SNPcat=='xMfam-xt',c(3:251)]<-
		sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:251)]-
		sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
	
	sampleFreq[SNPyr==year & SNPcat=='xFfam-xt',c(3:251)]<-
	  sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:251)]-
	  sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:251)]
}

#calculate variances and covariances
for(year in c(1999:2013)){
	sampleVar[varYr==year & varCat=='xMfam',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:251)])^2)
	
	sampleVar[varYr==year & varCat=='xFfam',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:251)])^2)
		
	sampleVar[varYr==year & varCat=='xMmend',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMmend',c(3:251)])^2)
	
	sampleVar[varYr==year & varCat=='xFmend',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFmend',c(3:251)])^2)
		
	sampleVar[varYr==year & varCat=='xMfam-xt',3]<-
		mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMfam-xt',c(3:251)])^2)
	
	sampleVar[varYr==year & varCat=='xFfam-xt',3]<-
	  mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFfam-xt',c(3:251)])^2)
}

#get date script is run
today<-format(Sys.Date(),format="%d%b%Y")

#save output
save(samplePars,sampleFreq,file=paste("modelZIntermediateFiles_",today,".rdata",sep=''))
save(sampleVar,file=paste("sampleVar_Z_SR",today,".rdata",sep=''))
