#script to estimate variance in variance in allele frequency change over time by bootstrapping  for Z loci
#Nancy Chen & Graham Coop & Rose Driscoll & Felix Beaudry
#Last updated: 31 May 2022

library(plyr)
library(tidyr)
library(foreach)
library(dplyr)
library(doParallel)
library(data.table)
`%notin%` <- Negate(`%in%`)
`%ni%` <- Negate(`%in%`)


####get input files####
load(file='../working_files/intermediate_files/indivlistgeno_Zpruned.rdata')
indivlistgeno <- indivlistgenoZ_LD[,-c(7)]
colnames(indivlistgeno)[1:7] <- c("Year","Indiv", "Category", "Genotyped", "Mom", "Dad", "Sex")

#Z SNP info
map <- fread('../working_files/geno.pruned.map')
ZSNPS <- map[map$V1 == 39 ,]
names(ZSNPS) <- c('Chr', 'SNP', 'genMap', 'SNPpos')
ZSNPS$SNP_name <- names(indivlistgeno[,-c(1:7)])

chip <- fread('../../genomeContent/genome_files/snp_annotation.txt', fill = TRUE, header = TRUE)
ZSNPS_chip <- left_join(ZSNPS, chip, by = c("SNP" = "SNP"))

#which scaffolds to keep
sizes <- fread('../working_files/FSJ.chrom.sizes')
Z_size <- sizes$V2[sizes$V1 %in% unique(ZSNPS_chip$scaffold)] 

w_size = 3400000
nloci<- 5000 #set sim loci number relative to window size; regular sim nloci / genome size in mb
#nloci=4

win_global = 0
for(win in seq(from=0,to=max(ZSNPS_chip$SNPpos,na.rm=T),by = w_size)){
  if(win+w_size<Z_size){
    #cat( ((win/w_size) + win_global)," ",win," ",Z_size,"\n")
    ZSNPS_chip$bootstrap[ZSNPS_chip$SNPpos > win & ZSNPS_chip$SNPpos <= (win+w_size) ] <- (win/w_size) + win_global
  }else{
    cat( ((win/w_size) + win_global)," ",win," ",Z_size," skipped\n")
    
  }
}


markersInPedOrder_chip_Z <- ZSNPS_chip[,c('SNP_name','bootstrap')]
names(markersInPedOrder_chip_Z)[1]<- "SNP"

markersInPedOrder_chip_Z_tally <- markersInPedOrder_chip_Z %>% group_by(bootstrap) %>% tally()
big_boots <- markersInPedOrder_chip_Z_tally$bootstrap[markersInPedOrder_chip_Z_tally$n >= 5]
markersInPedOrder_chip_Z <- markersInPedOrder_chip_Z[markersInPedOrder_chip_Z$bootstrap %in% big_boots,]


###
#start sim file
simindivlist <- indivlistgeno[order(indivlistgeno$Year),c(1:7)]
colnames(simindivlist) <- c( "Year","Indiv", "Category", "Genotyped", "Mom", "Dad", "Sex")
simindivgeno<-simindivlist[!duplicated(simindivlist$Indiv),]
colnames(simindivgeno) <- c( "Year","Indiv", "Category", "Genotyped", "Mom", "Dad", "Sex")

#separate into moms vs dads vs nestlings
simindivgenoMoms<-
  simindivgeno[(simindivgeno$Category!='nestling' | is.na(simindivgeno$Mom)) & simindivgeno$Sex==2,]
simindivgenoDads<-
  simindivgeno[(simindivgeno$Category!='nestling' | is.na(simindivgeno$Mom)) & simindivgeno$Sex==1,]
simindivgenoNestlings<-
  simindivgeno[simindivgeno$Category=='nestling' & !is.na(simindivgeno$Mom),]


#start sample file
allVar <- 
  rbind.data.frame(
    data.frame(Year=rep(c(1999:2013),each=28),Category=rep(c(
      'xt1-xt',
      'xMs-xt','xFs-xt',
      'xMi-xt','xFi-xt',
      'xMb-xt','xFb-xt',
      'xMsxMi','xFsxFi','xMsxMb','xFsxFb','xMixMb','xFixFb',
      'xMsxFs','xMixFi','xMbxFb',
      'xMsxFi','xMsxFb','xMixFb',
      'xFsxMi','xFsxMb','xFixMb',
      'xMfam','xFfam','xMmend','xFmend','xMfam-xt','xFfam-xt'),15),stringsAsFactors=FALSE)
    ,
    data.frame(Year=rep(c(1999:2013),each=112),Category=rep(c(
      'pt1-pt',
      'pMs-pt','sxMs-xt','errMS-errT','pMspterrMSerrT',
      'pFs-pt','sxFs-xt','errFS-errT','pFspterrFSerrT',
      'pMi-pt','sxMi-xt','errMI-errT','pMipterrMIerrT',
      'pFi-pt','sxFi-xt','errFI-errT','pFipterrFIerrT',
      'pMb-pt','sxMb-xt','errMB-errT','pMbpterrMBerrT',
      'pFb-pt','sxFb-xt','errFB-errT','pFbpterrFBerrT',
      'pMspMi','sxMsxMi','sxMserrMI','errMSxMi','errMSerrMI',
      'pMspMb','sxMsxMb','sxMserrMB','errMSxMb','errMSerrMB',
      'pMipMb','sxMixMb','sxMierrMB','errMIxMb','errMIerrMB',
      'pFspFi','sxFsxFi','sxFserrFI','errFSxFi','errFSerrFI',
      'pFspFb','sxFsxFb','sxFserrFB','errFSxFb','errFSerrFB',
      'pFipFb','sxFixFb','sxFierrFB','errFIxFb','errFIerrFB',
      'pMspFs','sxMsxFs','sxMserrFS','errMSxFs','errMSerrFS',
      'pMipFi','sxMixFi','sxMierrFI','errMIxFi','errMIerrFI',
      'pMbpFb','sxMbxFb','sxMberrFB','errMBxFb','errMBerrFB',
      'pMspFi','sxMsxFi','sxMserrFI','errMSxFi','errMSerrFI',
      'pMspFb','sxMsxFb','sxMserrFB','errMSxFb','errMSerrFB',
      'pMipFb','sxMixFb','sxMierrFB','errMIxFb','errMIerrFB',
      'pFspMi','sxFsxMi','sxFserrMI','errFSxMi','errFSerrMI',
      'pFspMb','sxFsxMb','sxFserrMB','errFSxMb','errFSerrMB',
      'pFipMb','sxFixMb','sxFierrMB','errFIxMb','errFIerrMB',
      'pMmend','sxMmend','errMMEND','pMfam-pt','sxMfam-xt','errMFAM-errT',
      'pFmend','sxFmend','errFMEND','pFfam-pt','sxFfam-xt','errFFAM-errT')
      ,15),stringsAsFactors=FALSE)
  )


####start the loop####
loop=1
#win=4
for(win in unique(na.omit(markersInPedOrder_chip_Z$bootstrap))){
  cat("window = ",win,"\n")
  
  #n = number of Genotyped individuals, x = sample allele frequency
  #calculate allele freq
  sampleFreq<-data.frame(Year=rep(c(1999:2013),each=31),Category=rep(c(
    'nt','xt','xt1-xt',
    'nMs','xMs','nFs','xFs', 'xMs-xt','xFs-xt',
    'nMi','xMi','nFi','xFi', 'xMi-xt','xFi-xt',
    'nMb','xMb','nFb','xFb', 'xMb-xt','xFb-xt',
    'xMdad','xFdad','xMmom','xFmom','xMfam','xFfam','xMmend','xFmend','xMfam-xt','xFfam-xt'),15),
    stringsAsFactors=FALSE)
  # 'xMdad', 'xFdad', 'xMmom', and 'xFmom' are the frequencies with which males transmit an 
  # allele to their sons, males transmit an allele to their daughters, females transmit an 
  # allele to their sons, and females transmit an allele to their daughters respectively
  # Omit xFmom for the Z as females never transmit their Z to their daughters
  
  sampleFreq<-rbind(data.frame(Year=rep(1998,2),Category=c('nt','xt'),
                               stringsAsFactors=FALSE),sampleFreq)
  SNPyr<-sampleFreq$Year
  SNPcat<-sampleFreq$Category
  
  igYear<-indivlistgeno$Year
  
  markers <- markersInPedOrder_chip_Z$SNP[markersInPedOrder_chip_Z$bootstrap == win]
  markers <- c(na.omit(markers))
    
  for(snp in markers){ 
    sampleFreq[SNPyr==1998 & SNPcat=='nt',snp]<-
      2*sum(!is.na(indivlistgeno[igYear==1998&indivlistgeno$Sex==1,snp])) +
      sum(!is.na(indivlistgeno[igYear==1998&indivlistgeno$Sex==2,snp]))
    
    sampleFreq[SNPyr==1998 & SNPcat=='xt',snp]<-
      sum(indivlistgeno[igYear==1998,snp],na.rm=TRUE)/
      sampleFreq[SNPyr==1998 & SNPcat=='nt',snp]
    
    for(year in c(1999:2013)){
      genoYr<-indivlistgeno[igYear==year,]
      
      sampleFreq[SNPyr==year & SNPcat=='nt',snp]<-
        2*sum(!is.na(genoYr[genoYr$Sex==1,snp])) + sum(!is.na(genoYr[genoYr$Sex==2,snp]))
      
      sampleFreq[SNPyr==year & SNPcat=='xt',snp]<-sum(genoYr[,snp],na.rm=TRUE)/
        sampleFreq[SNPyr==year & SNPcat=='nt',snp]
      
      sampleFreq[SNPyr==year & SNPcat=='nMs',snp]<-
        2*sum(!is.na(genoYr[genoYr$Category=='survivor'&genoYr$Sex==1,snp]))
      sampleFreq[SNPyr==year & SNPcat=='nFs',snp]<-
        sum(!is.na(genoYr[genoYr$Category=='survivor'&genoYr$Sex==2,snp]))
      
      sampleFreq[SNPyr==year & SNPcat=='xMs',snp]<-
        sum(genoYr[genoYr$Category=='survivor'&genoYr$Sex==1,snp],na.rm=TRUE)/
        sampleFreq[SNPyr==year & SNPcat=='nMs',snp]
      sampleFreq[SNPyr==year & SNPcat=='xFs',snp]<-
        sum(genoYr[genoYr$Category=='survivor'&genoYr$Sex==2,snp],na.rm=TRUE)/
        sampleFreq[SNPyr==year & SNPcat=='nFs',snp]
      
      sampleFreq[SNPyr==year & SNPcat=='nMi',snp]<-
        2*sum(!is.na(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==1,snp]))
      sampleFreq[SNPyr==year & SNPcat=='nFi',snp]<-
        sum(!is.na(genoYr[genoYr$Category=='immigrant'&genoYr$Sex==2,snp]))
      
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
        sum(!is.na(genoYr[genoYr$Category=='nestling'&genoYr$Sex==2,snp]))
      
      sampleFreq[SNPyr==year & SNPcat=='xMb',snp]<-
        sum(genoYr[genoYr$Category=='nestling'&genoYr$Sex==1,snp],na.rm=TRUE)/
        sampleFreq[SNPyr==year & SNPcat=='nMb',snp]
      sampleFreq[SNPyr==year & SNPcat=='xFb',snp]<-
        sum(genoYr[genoYr$Category=='nestling'&genoYr$Sex==2,snp],na.rm=TRUE)/
        sampleFreq[SNPyr==year & SNPcat=='nFb',snp]
    }
  }
  
  #calculate allele freq differences between each Category and the year before
  for(year in c(1999:2013)){
    #allele freq differences
    
    sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xt',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
    sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xt',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
    
    sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xMs',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
    sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xFs',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
    
    sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xMi',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
    sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xFi',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
    
    sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xMb',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
    sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xFb',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
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
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xt1-xt',c(3:length(sampleFreq))])^2)
    
    sampleVar[varYr==year & varCat=='xMs-xt','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:length(sampleFreq))])^2)
    sampleVar[varYr==year & varCat=='xFs-xt','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:length(sampleFreq))])^2)
    
    sampleVar[varYr==year & varCat=='xMi-xt','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:length(sampleFreq))])^2)
    sampleVar[varYr==year & varCat=='xFi-xt','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:length(sampleFreq))])^2)
    
    sampleVar[varYr==year & varCat=='xMb-xt','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:length(sampleFreq))])^2)
    sampleVar[varYr==year & varCat=='xFb-xt','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:length(sampleFreq))])^2)
    
    #covar
    sampleVar[varYr==year & varCat=='xMsxMi','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xFsxFi','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:length(sampleFreq))]))
    
    sampleVar[varYr==year & varCat=='xMsxMb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xFsxFb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:length(sampleFreq))]))
    
    sampleVar[varYr==year & varCat=='xMixMb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xFixFb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:length(sampleFreq))]))
    
    sampleVar[varYr==year & varCat=='xMsxFs','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xMixFi','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xMbxFb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:length(sampleFreq))]))
    
    sampleVar[varYr==year & varCat=='xMsxFi','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xMsxFb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xMixFb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFb-xt',c(3:length(sampleFreq))]))
    
    
    sampleVar[varYr==year & varCat=='xFsxMi','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMi-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xFsxMb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFs-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:length(sampleFreq))]))
    sampleVar[varYr==year & varCat=='xFixMb','avg']<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFi-xt',c(3:length(sampleFreq))])*
             as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMb-xt',c(3:length(sampleFreq))]))
    
  }
  
  #now calculate Mendelian noise
  #get unique individuals
  genoUnique<-indivlistgeno[!duplicated(indivlistgeno$Indiv),]
  
  #get sample allele frequencies of parents
  # year =1999
  for(year in c(1999:2013)){
    dadsofmales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==1,'Dad']
    dadsofmales<-data.frame(Indiv=dadsofmales[!is.na(dadsofmales)],stringsAsFactors=FALSE)
    dadsofmalesgeno<-merge(dadsofmales,genoUnique[,c(2,8:length(genoUnique))],by='Indiv',all.x=TRUE)
    
    momsofmales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==1,'Mom']
    momsofmales<-data.frame(Indiv=momsofmales[!is.na(momsofmales)],stringsAsFactors=FALSE)
    momsofmalesgeno<-merge(momsofmales,genoUnique[,c(2,8:length(genoUnique))],by='Indiv',all.x=TRUE)
    
    dadsoffemales<-indivlistgeno[indivlistgeno$Year==year & indivlistgeno$Category=='nestling' & indivlistgeno$Sex==2,'Dad']
    dadsoffemales<-data.frame(Indiv=dadsoffemales[!is.na(dadsoffemales)],stringsAsFactors=FALSE)
    dadsoffemalesgeno<-merge(dadsoffemales,genoUnique[,c(2,8:length(genoUnique))],by='Indiv',all.x=TRUE)
    
    #Moms don't contribute to daughters on Z so skip these
    #momsoffemales<-indivlist[indivlist$Year==year & indivlist$Category=='nestling' & indivlist$Sex==2,'Mom']
    #momsoffemales<-data.frame(Indiv=momsoffemales[!is.na(momsoffemales)],stringsAsFactors=FALSE)
    #momsoffemalesgeno<-merge(momsoffemales,genoUnique[,c(1,8:256)],by='Indiv',all.x=TRUE)
    
    for(snp in markers)
    {
      sampleFreq[SNPyr==year & SNPcat=='xMdad',snp]<-mean(dadsofmalesgeno[,snp],na.rm=TRUE)/2
      sampleFreq[SNPyr==year & SNPcat=='xMmom',snp]<-mean(momsofmalesgeno[,snp],na.rm=TRUE)
      
      sampleFreq[SNPyr==year & SNPcat=='xFdad',snp]<-mean(dadsoffemalesgeno[,snp],na.rm=TRUE)/2
      #sampleFreq[SNPyr==year & SNPcat=='xFmom',snp]<-mean(momsoffemalesgeno[,snp],na.rm=TRUE)/2
    }
  }
  
  #calculate other terms
  for(year in c(1999:2013)){
    sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:length(sampleFreq))]<-0.5*
      (sampleFreq[SNPyr==year & SNPcat=='xMdad',c(3:length(sampleFreq))]+
         sampleFreq[SNPyr==year & SNPcat=='xMmom',c(3:length(sampleFreq))]) 
    #here we multiply mom by 0.5 since her contribution makes up only 1/2 of her son's genotype
    
    sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xFdad',c(3:length(sampleFreq))]
    #here we don't multiply dad by 0.5 since his contribution makes up 100% of his daughter's genotype
    
    sampleFreq[SNPyr==year & SNPcat=='xMmend',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xMb',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:length(sampleFreq))]
    
    sampleFreq[SNPyr==year & SNPcat=='xFmend',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xFb',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:length(sampleFreq))]
    
    sampleFreq[SNPyr==year & SNPcat=='xMfam-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
    
    sampleFreq[SNPyr==year & SNPcat=='xFfam-xt',c(3:length(sampleFreq))]<-
      sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:length(sampleFreq))]-
      sampleFreq[SNPyr==(year-1) & SNPcat=='xt',c(3:length(sampleFreq))]
  }
  
  #calculate variances and covariances
  for(year in c(1999:2013)){
    sampleVar[varYr==year & varCat=='xMfam',3]<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMfam',c(3:length(sampleFreq))])^2)
    
    sampleVar[varYr==year & varCat=='xFfam',3]<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFfam',c(3:length(sampleFreq))])^2)
    
    sampleVar[varYr==year & varCat=='xMmend',3]<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMmend',c(3:length(sampleFreq))])^2)
    
    sampleVar[varYr==year & varCat=='xFmend',3]<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFmend',c(3:length(sampleFreq))])^2)
    
    sampleVar[varYr==year & varCat=='xMfam-xt',3]<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xMfam-xt',c(3:length(sampleFreq))])^2)
    
    sampleVar[varYr==year & varCat=='xFfam-xt',3]<-
      mean(as.numeric(sampleFreq[SNPyr==year & SNPcat=='xFfam-xt',c(3:length(sampleFreq))])^2)
  }
  ####start sim####
  
  #get real frequency of each allele in 1990 (accounting for different total # of alleles in males & females)
  datafreq1990<-laply(markers,function(x) 
    sum(indivlistgeno[indivlistgeno$Year==1990,x],na.rm=TRUE)/
      ((2*sum(!is.na(indivlistgeno[indivlistgeno$Year==1990&indivlistgeno$Sex==1,x])))
       +(sum(!is.na(indivlistgeno[indivlistgeno$Year==1990&indivlistgeno$Sex==2,x])))))
  
  #randomly sample from real allele frequencies
  simfreq<-sample(datafreq1990,nloci,replace=TRUE)
  
  
  #simulate genotypes for adults
  #moms - 0 1
  num.moms<-nrow(simindivgenoMoms)
  mom.genos<-sapply(1:nloci,function(loc){
    freq<-simfreq[loc]
    HWE<-c(1-freq,freq)	
    loc.genos<-sample(0:1,size=num.moms,prob=HWE,replace=TRUE)
    loc.genos
  })
  #dads - 0 1 2
  num.dads<-nrow(simindivgenoDads)
  dad.genos<-sapply(1:nloci,function(loc){
    freq<-simfreq[loc]
    HWE<-c((1-freq)^2,2*freq*(1-freq),freq^2)	
    loc.genos<-sample(0:2,size=num.dads,prob=HWE,replace=TRUE)
    loc.genos
  })
  
  #combine simulated mom genotypes with real mom IDs & info
  simindivgenoMoms<-cbind(simindivgenoMoms,mom.genos)
  #combine simulated dad genotypes with real dad IDs & info
  simindivgenoDads<-cbind(simindivgenoDads,dad.genos)
  
  #put the moms and dads back together in one table
  simindivgenoParents <- rbind(simindivgenoMoms, simindivgenoDads)
  
  #combine it back with a big table that has space for every individual
  simindivgenoAll<-cbind(simindivgeno,simindivgenoParents[match(simindivgeno$Indiv,
                                                                simindivgenoParents$Indiv),8:(nloci+7)])
  #add rownames
  rownames(simindivgenoAll)<-simindivgenoAll$Indiv
  
  #simulate genotypes for nestlings via Mendelian transmission of alleles from parents
  
  #get years
  nest.years<-unique(simindivgenoNestlings$Year)
  
  #for fathers, will need to sample their gametes to get kid genotypes
  make.male.gametes<-function(g){
    gametes<-rep(NA,length(g))
    gametes[g==0]<-0
    gametes[g==2]<-1
    hets <- g==1 & !is.na(g)
    gametes[hets]<-sample(c(0,1),sum(hets),prob=c(0.5,0.5),replace=TRUE)
    gametes
  }
  
  #get kid genotypes, one year at a time
  #year= 1991
  for(year in nest.years){
    #pull out female nestlings for this year
    these.female.nestlings<-simindivgenoNestlings[simindivgenoNestlings$Year==year 
                                                  & simindivgenoNestlings$Sex==2, "Indiv"]
    #pull out male nestlings for this year
    these.male.nestlings<-simindivgenoNestlings[simindivgenoNestlings$Year==year 
                                                & simindivgenoNestlings$Sex==1, "Indiv"]
    
    #ignore moms of this year's female nestlings as they don't contribute Zs
    #these.moms.of.daughters<-simindivgenoNestlings[simindivgenoNestlings$Year==year
    #& simindivgenoNestlings$Sex==2,'Mom']
    #these.moms.of.sons<-simindivgenoNestlings[simindivgenoNestlings$Year==year,'Mom']
    #these.moms.of.daughters<-simindivgenoNestlings$Mom[simindivgenoNestlings$Year==year
    #                                                   & simindivgenoNestlings$Sex==2]
    
    #get a list of the moms of this year's male nestlings
    these.moms.of.sons<-simindivgenoNestlings[simindivgenoNestlings$Year==year
                                              & simindivgenoNestlings$Sex==1, "Mom"]
    #get a list of the dads of this year's female nestlings
    these.dads.of.daughters<-simindivgenoNestlings[simindivgenoNestlings$Year==year
                                                   & simindivgenoNestlings$Sex==2, "Dad"]
    #get a list of the dads of this year's male nestlings
    these.dads.of.sons<-simindivgenoNestlings[simindivgenoNestlings$Year==year
                                              & simindivgenoNestlings$Sex==1, "Dad"]
    
    #check that the nestlings don't have genotypes yet (should be NA)
    stopifnot(all(is.na(simindivgenoAll[these.female.nestlings,8:(nloci+7)])))
    stopifnot(all(is.na(simindivgenoAll[these.male.nestlings,8:(nloci+7)])))
    
    #get parents' genotypes
    #moms.of.daughters.geno<-simindivgenoAll[these.moms.of.daughters,]
    moms.of.sons.geno<-simindivgenoAll[these.moms.of.sons,]
    dads.of.daughters.geno<-simindivgenoAll[these.dads.of.daughters,]
    dads.of.sons.geno<-simindivgenoAll[these.dads.of.sons,]
    
    #check that all parents are Genotyped (NOT NA)
    #stopifnot(all(!is.na(moms.of.daughters.geno[,8:(nloci+7)])))
    stopifnot(all(!is.na(moms.of.sons.geno[,8:(nloci+7)])))
    stopifnot(all(!is.na(dads.of.daughters.geno[,8:(nloci+7)])))
    stopifnot(all(!is.na(dads.of.sons.geno[,8:(nloci+7)])))
    
    #run the gamete selector function to pick which gamete dads give their kids
    Dads.of.daughters.gamete<-apply(dads.of.daughters.geno[,8:(nloci+7)],2,make.male.gametes)
    Dads.of.sons.gamete<-apply(dads.of.sons.geno[,8:(nloci+7)],2,make.male.gametes)
    # moms don't give daughters Z alleles so skip moms of daughters
    #get mom's genotype -> this is the gamete she gives her son
    Moms.of.sons.gamete<-moms.of.sons.geno[,8:(nloci+7)]
    
    #female nestlings just get the allele from their dad as their geno
    female.nestling.geno<-Dads.of.daughters.gamete
    #male nestlings get allele from dad + allele from mom
    male.nestling.geno<-Dads.of.sons.gamete+Moms.of.sons.gamete
    
    #stick the genotypes for this year's nestlings into the simindivgenoAll table
    simindivgenoAll[simindivgenoAll$Indiv %in% these.female.nestlings,8:(nloci+7)] <- female.nestling.geno
    simindivgenoAll[simindivgenoAll$Indiv %in% these.male.nestlings,8:(nloci+7)] <- male.nestling.geno
    
    #check that the nestlings now have genotypes (not NA)
    stopifnot(all(!is.na(simindivgenoAll[these.female.nestlings,8:(nloci+7)])))
    stopifnot(all(!is.na(simindivgenoAll[these.male.nestlings,8:(nloci+7)])))
    
  }
  
  #now we want to have each indiv appear multiple times again
  simdataTrue<-merge(simindivlist,simindivgenoAll[,c(2,8:(nloci+7))], by.x='Indiv',by.y='Indiv',all.x=TRUE)	
  
  
  
  
  #calculate sample allele freq
  #mimic sampling of Genotyped indiv by selecting only indivs who actually were Genotyped
  simdataSample<-simdataTrue[simdataTrue$Genotyped=='Y',]
  #simdataSample[c(1:20),c(1:20)]
  #get unique indivs in simulated data (all & Genotyped)
  simdataTrueUnique<-simdataTrue[!duplicated(simdataTrue$Indiv),]
  simdataSampleUnique<-simdataSample[!duplicated(simdataSample$Indiv),]
  
  #calculate population (p) and sample allele freq (x)
  #and the error in allele freq estimation due to sampling: err = x-p
  
  #create data frame to hold simulated allele freqs
  simAlleleFreq<-data.frame(Year=integer(),Category=character(),stringsAsFactors=FALSE)
  
  year<-1998
  sim<-foreach(i=names(simdataTrue)[8:(nloci+7)],.combine=cbind) %do% {
    #create data frame to hold allele freqs & error
    tmp<-data.frame(Year=rep(year,each=3),Category=c('pt','sxt','errT'),
                    stringsAsFactors=FALSE)
    
    frqYr1<-tmp$Year
    frqCat1<-tmp$Category
    
    tmp[frqYr1==year & frqCat1=='pt',3]<-
      sum(simdataTrue[simdataTrue$Year==year,i],na.rm=TRUE)/
      ((2*sum(!is.na(simdataTrue[simdataTrue$Year==year&simdataTrue$Sex==1,i])))
       +(sum(!is.na(simdataTrue[simdataTrue$Year==year&simdataTrue$Sex==2,i]))))
    
    tmp[frqYr1==year & frqCat1=='sxt',3]<-
      sum(simdataSample[simdataSample$Year==year,i],na.rm=TRUE)/
      ((2*sum(!is.na(simdataSample[simdataSample$Year==year&simdataSample$Sex==1,i])))
       +(sum(!is.na(simdataSample[simdataSample$Year==year&simdataSample$Sex==2,i]))))
    
    tmp[,3]
  }
  
  
  #Names for the values we just calculated (year and Category/parameter)
  simName<-data.frame(Year=rep(year,each=3),Category=c('pt','sxt','errT'),
                      stringsAsFactors=FALSE)
  
  #Add the simulation results to the data frame
  #err rows are NA b/c we have not calculated the error yet
  sim1<-cbind(simName,sim)	
  #Add the simulation data to simAlleleFreq (we'll collect the data from all years here)
  simAlleleFreq<-rbind(simAlleleFreq,sim1)
  
  
  for(year in c(1999:2013)){
    #get moms of sons, dads of sons, and dads of daughters for this year
    
    #get moms of male nestlings born this year
    moms_of_sons<-simdataTrue[simdataTrue$Year==year & simdataTrue$Category=='nestling' & simdataTrue$Sex==1,'Mom']
    
    
    #convert list of moms of sons to a data frame
    moms_of_sons<-data.frame(Indiv=moms_of_sons[!is.na(moms_of_sons)],stringsAsFactors=FALSE)
    #collect simulated moms of sons genotypes (including those simulated for unGenotyped indivs) from simdataTrueUnique
    moms_of_sons_geno<-merge(moms_of_sons,simdataTrueUnique[,c(1,8:(nloci+7))],by.x='Indiv',by.y='Indiv',
                             all.x=TRUE)
    #collect simulated moms of sons genotypes (sampled based on real genotyping status) from simdataSampleUnique
    moms_of_sons_genoSample<-merge(moms_of_sons,simdataSampleUnique[,c(1,8:(nloci+7))],by.x='Indiv',
                                   by.y='Indiv',all.x=TRUE)
    #many of the rows in this table are NAs because of the sampling
    
    #don't need moms of daughters as moms do not contribute to daughters
    
    #get dads of male nestlings born this year
    dads_of_sons<-simdataTrue[simdataTrue$Year==year & simdataTrue$Category=='nestling' & simdataTrue$Sex==1,'Dad']
    #convert list of dads of sons to a data frame
    dads_of_sons<-data.frame(Indiv=dads_of_sons[!is.na(dads_of_sons)],stringsAsFactors=FALSE)
    #collect simulated dads of sons genotypes (including those simulated for unGenotyped indivs) from simdataTrueUnique
    dads_of_sons_geno<-merge(dads_of_sons,simdataTrueUnique[,c(1,8:(nloci+7))],by.x='Indiv',by.y='Indiv',
                             all.x=TRUE)
    #collect simulated dads of sons genotypes (sampled based on real genotyping status) from simdataSampleUnique
    dads_of_sons_genoSample<-merge(dads_of_sons,simdataSampleUnique[,c(1,8:(nloci+7))],by.x='Indiv',
                                   by.y='Indiv',all.x=TRUE)
    #many of the rows in this table are NAs because of the sampling
    
    #get dads of female nestlings born this year
    dads_of_daughters<-simdataTrue[simdataTrue$Year==year & simdataTrue$Category=='nestling' & simdataTrue$Sex==2,'Dad']
    #convert list of dads of daughters to a data frame
    dads_of_daughters<-data.frame(Indiv=dads_of_daughters[!is.na(dads_of_daughters)],stringsAsFactors=FALSE)
    #collect simulated dads of daughters genotypes (including those simulated for unGenotyped indivs) from simdataTrueUnique
    dads_of_daughters_geno<-merge(dads_of_daughters,simdataTrueUnique[,c(1,8:(nloci+7))],by.x='Indiv',by.y='Indiv',
                                  all.x=TRUE)
    #collect simulated dads of daughters genotypes (sampled based on real genotyping status) from simdataSampleUnique
    dads_of_daughters_genoSample<-merge(dads_of_daughters,simdataSampleUnique[,c(1,8:(nloci+7))],by.x='Indiv',
                                        by.y='Indiv',all.x=TRUE)
    #many of the rows in this table are NAs because of the sampling
    
    #for each snp
    sim<-foreach(i=names(simdataTrue)[8:(nloci+7)],.combine=cbind) %do% {
      #make a data frame to put all these parameters in for each year
      tmp<-data.frame(Year=rep(year,each=70),Category=c(
        'pt','sxt','errT', 'pt1-pt',
        
        'pMs','sxMs','errMS', 'pMs-pt', 'sxMs-xt', 'errMS-errT',
        'pMi','sxMi','errMI', 'pMi-pt','sxMi-xt','errMI-errT',
        'pMb','sxMb','errMB','pMb-pt','sxMb-xt','errMB-errT',
        'pMdad','sxMdad','errMdad', 
        'pMmom','sxMmom','errMmom',
        'pMfam','sxMfam','pMfam-pt', 'sxMfam-xt', 'errMFAM', 'errMFAM-errT',
        'pMmend','sxMmend','errMMEND', 
        
        'pFs','sxFs','errFS', 'pFs-pt','sxFs-xt','errFS-errT',
        'pFi','sxFi','errFI', 'pFi-pt','sxFi-xt','errFI-errT',
        'pFb','sxFb','errFB', 'pFb-pt','sxFb-xt','errFB-errT',
        'pFdad','sxFdad','errFdad',
        'pFmom','sxFmom','errFmom',
        'pFfam','sxFfam','pFfam-pt','sxFfam-xt', 'errFFAM',  'errFFAM-errT',
        'pFmend','sxFmend','errFMEND'
      ),stringsAsFactors=FALSE)
      #previously pm = allele frequency of all mothers, pf = allele frequency of all fathers,
      #and similar for xm, xf, errM, errF
      #Now instead we are using pMdad = allele frequency of all fathers of sons,
      #pMmom = allele frequency of all mothers of sons,
      #pFdad = allele frequency of all fathers of daughters,
      #and likewise xMdad, xMmom, errMdad, errMmom, etc.
      
      
      frqYr1<-tmp$Year
      frqCat1<-tmp$Category
      
      #pt is just the mean of all of the simulated data (no sampling)
      tmp[frqYr1==year & frqCat1=='pt',3]<-
        sum(simdataTrue[simdataTrue$Year==year,i],na.rm=TRUE)/
        ((2*sum(!is.na(simdataTrue[simdataTrue$Year==year&simdataTrue$Sex==1,i])))
         +(sum(!is.na(simdataTrue[simdataTrue$Year==year&simdataTrue$Sex==2,i]))))
      
      #xt is the mean of the sampled simulated data
      tmp[frqYr1==year & frqCat1=='sxt',3]<-
        sum(simdataSample[simdataSample$Year==year,i],na.rm=TRUE)/
        ((2*sum(!is.na(simdataSample[simdataSample$Year==year&simdataSample$Sex==1,i])))
         +(sum(!is.na(simdataSample[simdataSample$Year==year&simdataSample$Sex==2,i]))))
      
      #ps
      tmp[frqYr1==year & frqCat1=='pMs',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                               simdataTrue$Category=='survivor' & simdataTrue$Sex==1,i])/2
      tmp[frqYr1==year & frqCat1=='pFs',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                               simdataTrue$Category=='survivor' & simdataTrue$Sex==2,i])
      
      #xs
      tmp[frqYr1==year & frqCat1=='sxMs',3]<-mean(simdataSample[simdataSample$Year==year & 
                                                                  simdataSample$Category=='survivor' & simdataSample$Sex==1,i])/2
      tmp[frqYr1==year & frqCat1=='sxFs',3]<-mean(simdataSample[simdataSample$Year==year & 
                                                                  simdataSample$Category=='survivor' & simdataSample$Sex==2,i])
      
      #pi
      tmp[frqYr1==year & frqCat1=='pMi',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                               simdataTrue$Category=='immigrant' & simdataTrue$Sex==1,i])/2
      tmp[frqYr1==year & frqCat1=='pFi',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                               simdataTrue$Category=='immigrant' & simdataTrue$Sex==2,i])
      
      #xi
      tmp[frqYr1==year & frqCat1=='sxMi',3]<-
        ifelse(is.na(mean(simdataSample[simdataSample$Year==year & 
                                          simdataSample$Category=='immigrant' & simdataSample$Sex==1,i])),0,
               mean(simdataSample[simdataSample$Year==year & 
                                    simdataSample$Category=='immigrant' & simdataSample$Sex==1,i])/2)
      tmp[frqYr1==year & frqCat1=='sxFi',3]<-
        ifelse(is.na(mean(simdataSample[simdataSample$Year==year & 
                                          simdataSample$Category=='immigrant' & simdataSample$Sex==2,i])),0,
               mean(simdataSample[simdataSample$Year==year & 
                                    simdataSample$Category=='immigrant' & simdataSample$Sex==2,i]))
      #ifelse here to catch years with no Genotyped imms
      #imms are the only Category that sometimes is 0 - we checked
      
      #pb
      tmp[frqYr1==year & frqCat1=='pMb',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                               simdataTrue$Category=='nestling' & simdataTrue$Sex==1,i])/2
      tmp[frqYr1==year & frqCat1=='pFb',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                               simdataTrue$Category=='nestling' & simdataTrue$Sex==2,i])
      
      #xb
      tmp[frqYr1==year & frqCat1=='sxMb',3]<-mean(simdataSample[simdataSample$Year==year &
                                                                  simdataSample$Category=='nestling' & simdataSample$Sex==1,i])/2
      tmp[frqYr1==year & frqCat1=='sxFb',3]<-mean(simdataSample[simdataSample$Year==year &
                                                                  simdataSample$Category=='nestling' & simdataSample$Sex==2,i])
      
      #pMmom & xMmom
      tmp[frqYr1==year & frqCat1=='pMmom',3]<-mean(moms_of_sons_geno[,i],na.rm=TRUE)
      tmp[frqYr1==year & frqCat1=='sxMmom',3]<-mean(moms_of_sons_genoSample[,i],na.rm=TRUE)
      
      #pMdad & xMdad
      tmp[frqYr1==year & frqCat1=='pMdad',3]<-mean(dads_of_sons_geno[,i],na.rm=TRUE)/2
      tmp[frqYr1==year & frqCat1=='sxMdad',3]<-mean(dads_of_sons_genoSample[,i],na.rm=TRUE)/2
      
      #pFdad & xFdad
      tmp[frqYr1==year & frqCat1=='pFdad',3]<-mean(dads_of_daughters_geno[,i],na.rm=TRUE)/2
      tmp[frqYr1==year & frqCat1=='sxFdad',3]<-mean(dads_of_daughters_genoSample[,i],na.rm=TRUE)/2
      
      tmp[,3]
    }
    
    #categories (parameter names) and years to combine with the results of our calculations
    simName<-data.frame(Year=rep(year,each=70),Category=c(
      'pt','sxt','errT', 'pt1-pt',
      
      'pMs','sxMs','errMS', 'pMs-pt', 'sxMs-xt', 'errMS-errT',
      'pMi','sxMi','errMI', 'pMi-pt','sxMi-xt','errMI-errT',
      'pMb','sxMb','errMB','pMb-pt','sxMb-xt','errMB-errT',
      'pMdad','sxMdad','errMdad', 
      'pMmom','sxMmom','errMmom',
      'pMfam','sxMfam','pMfam-pt', 'sxMfam-xt', 'errMFAM', 'errMFAM-errT',
      'pMmend','sxMmend','errMMEND', 
      
      'pFs','sxFs','errFS', 'pFs-pt','sxFs-xt','errFS-errT',
      'pFi','sxFi','errFI', 'pFi-pt','sxFi-xt','errFI-errT',
      'pFb','sxFb','errFB', 'pFb-pt','sxFb-xt','errFB-errT',
      'pFdad','sxFdad','errFdad',
      'pFmom','sxFmom','errFmom',
      'pFfam','sxFfam','pFfam-pt','sxFfam-xt', 'errFFAM',  'errFFAM-errT',
      'pFmend','sxFmend','errFMEND'
    ),stringsAsFactors=FALSE)
    #add the names to the calculation results
    sim1<-cbind(simName,sim)	
    
    #combine with simAlleleFreq which is where we are collecting the results from all the snps
    simAlleleFreq<-rbind(simAlleleFreq,sim1)
  }
  
  
  #calculate error and allele freq differences between each Category and the year before
  #err = true error (no hypergeometric error)
  frqYr<-simAlleleFreq$Year
  frqCat<-simAlleleFreq$Category
  
  #err = true error
  simAlleleFreq[frqYr==1998 & frqCat=='errT',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==1998 & frqCat=='sxt',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==1998 & frqCat=='pt',c(3:(nloci+2))]
  
  #calcuate error for each year based on difference b/w p and x
  for(year in c(1999:2013)){
    #true error
    simAlleleFreq[frqYr==year & frqCat=='errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxt',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pt',c(3:(nloci+2))]	
    
    simAlleleFreq[frqYr==year & frqCat=='errMS',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMs',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pMs',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='errFS',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFs',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pFs',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='errMI',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMi',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pMi',c(3:(nloci+2))]	
    simAlleleFreq[frqYr==year & frqCat=='errFI',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFi',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pFi',c(3:(nloci+2))]	
    
    simAlleleFreq[frqYr==year & frqCat=='errMB',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pMb',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='errFB',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pFb',c(3:(nloci+2))]
    
    #mom
    simAlleleFreq[frqYr==year & frqCat=='errMmom',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMmom',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pMmom',c(3:(nloci+2))]
    
    #dad
    simAlleleFreq[frqYr==year & frqCat=='errMdad',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMdad',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pMdad',c(3:(nloci+2))]	
    simAlleleFreq[frqYr==year & frqCat=='errFdad',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFdad',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pFdad',c(3:(nloci+2))]	
    
    
    #allele freq differences
    #true (population) freq (p)
    #this is actually doing pt - pt-1 even thought the Category name is pt1-pt
    #so the math does match my equations
    simAlleleFreq[frqYr==year & frqCat=='pt1-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pt',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pMs',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pFs',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pMi',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pFi',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pMb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pFb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    
    #sample freq (x)
    simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMs',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='sxt',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFs',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='sxt',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMi',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='sxt',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFi',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='sxt',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='sxt',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='sxt',c(3:(nloci+2))]
    
    #true errors
    simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errMS',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errFS',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errMI',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errFI',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errMB',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errFB',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
    
    #Mendelian noise
    simAlleleFreq[frqYr==year & frqCat=='pMfam',c(3:(nloci+2))]<-
      0.5*(simAlleleFreq[frqYr==year & frqCat=='pMmom',c(3:(nloci+2))]+
             simAlleleFreq[frqYr==year & frqCat=='pMdad',c(3:(nloci+2))])
    #here we also multiply mom by 0.5 since her contribution makes up only 1/2 of her son's genotype
    
    simAlleleFreq[frqYr==year & frqCat=='pFfam',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pFdad',c(3:(nloci+2))]
    #here we don't multiply dad by 0.5 since his contribution makes up 100% of his daughter's genotype
    
    
    simAlleleFreq[frqYr==year & frqCat=='pMmend',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pMb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pMfam',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='pFmend',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pFb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='pFfam',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='pMfam-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pMfam',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='pFfam-pt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='pFfam',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='sxMfam',c(3:(nloci+2))]<-
      0.5*(simAlleleFreq[frqYr==year & frqCat=='sxMmom',c(3:(nloci+2))]+
             simAlleleFreq[frqYr==year & frqCat=='sxMdad',c(3:(nloci+2))])
    #here we also multiply mom by 0.5 since her contribution makes up only 1/2 of her son's genotype
    
    simAlleleFreq[frqYr==year & frqCat=='sxFfam',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFdad',c(3:(nloci+2))]
    #here we don't multiply dad by 0.5 since his contribution makes up 100% of his daughter's genotype
    
    simAlleleFreq[frqYr==year & frqCat=='sxMmend',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='sxMfam',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='sxFmend',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFb',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='sxFfam',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='sxMfam-xt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxMfam',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='sxt',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='sxFfam-xt',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='sxFfam',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='sxt',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='errMFAM',c(3:(nloci+2))]<-
      0.5*(simAlleleFreq[frqYr==year & frqCat=='errMmom',c(3:(nloci+2))]+
             simAlleleFreq[frqYr==year & frqCat=='errMdad',c(3:(nloci+2))])
    #here we also multiply mom by 0.5 since her contribution makes up only 1/2 of her son's genotype
    
    simAlleleFreq[frqYr==year & frqCat=='errFFAM',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errFdad',c(3:(nloci+2))]
    #here we don't multiply dad by 0.5 since his contribution makes up 100% of his daughter's genotype
    
    simAlleleFreq[frqYr==year & frqCat=='errMMEND',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errMB',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='errMFAM',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='errFMEND',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errFB',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==year & frqCat=='errFFAM',c(3:(nloci+2))]
    
    simAlleleFreq[frqYr==year & frqCat=='errMFAM-errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errMFAM',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
    simAlleleFreq[frqYr==year & frqCat=='errFFAM-errT',c(3:(nloci+2))]<-
      simAlleleFreq[frqYr==year & frqCat=='errFFAM',c(3:(nloci+2))]-
      simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
    
  }
  
  
  #calculate variances and covariances
  simVar<-data.frame(Year=rep(c(1999:2013),each=112),Category=rep(c(
    'pt1-pt',
    'pMs-pt','sxMs-xt','errMS-errT','pMspterrMSerrT',
    'pFs-pt','sxFs-xt','errFS-errT','pFspterrFSerrT',
    'pMi-pt','sxMi-xt','errMI-errT','pMipterrMIerrT',
    'pFi-pt','sxFi-xt','errFI-errT','pFipterrFIerrT',
    'pMb-pt','sxMb-xt','errMB-errT','pMbpterrMBerrT',
    'pFb-pt','sxFb-xt','errFB-errT','pFbpterrFBerrT',
    'pMspMi','sxMsxMi','sxMserrMI','errMSxMi','errMSerrMI',
    'pMspMb','sxMsxMb','sxMserrMB','errMSxMb','errMSerrMB',
    'pMipMb','sxMixMb','sxMierrMB','errMIxMb','errMIerrMB',
    'pFspFi','sxFsxFi','sxFserrFI','errFSxFi','errFSerrFI',
    'pFspFb','sxFsxFb','sxFserrFB','errFSxFb','errFSerrFB',
    'pFipFb','sxFixFb','sxFierrFB','errFIxFb','errFIerrFB',
    'pMspFs','sxMsxFs','sxMserrFS','errMSxFs','errMSerrFS',
    'pMipFi','sxMixFi','sxMierrFI','errMIxFi','errMIerrFI',
    'pMbpFb','sxMbxFb','sxMberrFB','errMBxFb','errMBerrFB',
    'pMspFi','sxMsxFi','sxMserrFI','errMSxFi','errMSerrFI',
    'pMspFb','sxMsxFb','sxMserrFB','errMSxFb','errMSerrFB',
    'pMipFb','sxMixFb','sxMierrFB','errMIxFb','errMIerrFB',
    'pFspMi','sxFsxMi','sxFserrMI','errFSxMi','errFSerrMI',
    'pFspMb','sxFsxMb','sxFserrMB','errFSxMb','errFSerrMB',
    'pFipMb','sxFixMb','sxFierrMB','errFIxMb','errFIerrMB',
    'pMmend','sxMmend','errMMEND','pMfam-pt','sxMfam-xt','errMFAM-errT',
    'pFmend','sxFmend','errFMEND','pFfam-pt','sxFfam-xt','errFFAM-errT')
    ,15),stringsAsFactors=FALSE)
  
  bsYr<-simVar$Year
  bsCat<-simVar$Category
  for(year in c(1999:2013)){
    #total variance, for males and females
    simVar[bsYr==year & bsCat=='pt1-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pt1-pt',c(3:(nloci+2))])^2)
    
    #variance for each Category, for males and females
    #survivors
    simVar[bsYr==year & bsCat=='pMs-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxMs-xt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errMS-errT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='pMspterrMSerrT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pFs-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxFs-xt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errFS-errT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='pFspterrFSerrT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))]))
    
    #imms
    simVar[bsYr==year & bsCat=='pMi-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxMi-xt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errMI-errT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='pMipterrMIerrT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pFi-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxFi-xt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',3:(nloci+2)])^2)
    simVar[bsYr==year & bsCat=='errFI-errT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='pFipterrFIerrT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
    
    #births
    simVar[bsYr==year & bsCat=='pMb-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxMb-xt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errMB-errT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='pMbpterrMBerrT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pFb-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxFb-xt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errFB-errT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='pFbpterrFBerrT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    
    
    #covariance between categories - males
    simVar[bsYr==year & bsCat=='pMspMi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMsxMi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMserrMI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSxMi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSerrMI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pMspMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMsxMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMserrMB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSxMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSerrMB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pMipMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMixMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMierrMB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMIxMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMIerrMB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    
    #covariance between categories - females
    simVar[bsYr==year & bsCat=='pFspFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFsxFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFserrFI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFSxFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFSerrFI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pFspFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFsxFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFserrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFSxFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFSerrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pFipFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFixFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFierrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFIxFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFIerrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    
    
    #covariance between males and females within a Category
    #survivors
    simVar[bsYr==year & bsCat=='pMspFs',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMsxFs',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMserrFS',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSxFs',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSerrFS',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))]))
    
    #imms
    simVar[bsYr==year & bsCat=='pMipFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMixFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMierrFI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMIxFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMIerrFI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
    
    #births
    simVar[bsYr==year & bsCat=='pMbpFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMbxFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMberrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMBxFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMBerrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    
    
    #covariance between males and females in different categories
    simVar[bsYr==year & bsCat=='pMspFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMsxFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMserrFI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSxFi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSerrFI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pMspFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMsxFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMserrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSxFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMSerrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pMipFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMixFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxMierrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMIxFb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errMIerrFB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pFspMi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFsxMi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFserrMI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFSxMi',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMi-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFSerrMI',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pFspMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFsxMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFserrMB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFs-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFSxMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFSerrMB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    
    simVar[bsYr==year & bsCat=='pFipMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFixMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='sxFierrMB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFi-xt',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFIxMb',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMb-xt',c(3:(nloci+2))]))
    simVar[bsYr==year & bsCat=='errFIerrMB',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])*
             as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
    
    
    #variance for mendelian segregation and family size 
    #males
    simVar[bsYr==year & bsCat=='pMmend',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMmend',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxMmend',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMmend',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errMMEND',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMMEND',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='pMfam-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMfam-pt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxMfam-xt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxMfam-xt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errMFAM-errT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMFAM-errT',c(3:(nloci+2))])^2)
    
    #females
    simVar[bsYr==year & bsCat=='pFmend',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFmend',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxFmend',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFmend',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errFMEND',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFMEND',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='pFfam-pt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFfam-pt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='sxFfam-xt',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='sxFfam-xt',c(3:(nloci+2))])^2)
    simVar[bsYr==year & bsCat=='errFFAM-errT',3]<-
      mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFFAM-errT',c(3:(nloci+2))])^2)
  }
  
  
  names(simVar) <- names(sampleVar)
  
  allVar_tmp <- rbind.data.frame(sampleVar,simVar)
  allVar[,(loop+2)] <- allVar_tmp[,3]
  names(allVar)[(loop+2)] <- paste("bs",win,sep="")
  loop = loop +1
  today<-format(Sys.Date(),format="%d%b%Y")
  save(allVar,file=paste("allVar_int_boot_Z_w",(w_size/1000000),"mb_",today,"_pruned.rdata",sep=''))
  
}

#save output
today<-format(Sys.Date(),format="%d%b%Y")
save(allVar,file=paste("allVar_boot_Z_w",(w_size/1000000),"mb_",today,"_pruned.rdata",sep=''))










