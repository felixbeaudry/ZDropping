#script to model variance in allele frequency change over time
#run sims for n SNPs to get empirical estimates of covariances and errors for Z loci
#Nancy Chen & Graham Coop & Rose Driscoll & Felix Beaudry
#Last updated: Jul 7 2021

library(plyr)
library(foreach)
library(doParallel)
library(data.table)
library(dplyr)

####set variables and make/import tables####
#number of SNPs to simulate
nloci<-100000
#nloci=2

#get date script is run
today<-format(Sys.Date(),format="%d%b%Y")

#get input files: fixed list of indiv in each Category each year
load(file='~/Downloads/genedropZplinkInput/indivlistgenoZ_26Apr2022.rdata')
indivlistgeno <- indivlistgenoZ


#indivlistgeno$Indiv<-as.character(indivlistgeno$Indiv)
#indivlistgeno$Dad<-as.character(indivlistgeno$Dad)
#indivlistgeno$Mom<-as.character(indivlistgeno$Mom)

####simulate starting genotypes####
#get real frequency of each allele in 1990 (accounting for different total # of alleles in males & females)
datafreq1990<-laply(names(indivlistgeno)[9:length(indivlistgeno)],function(x) 
  sum(indivlistgeno[indivlistgeno$Year==1990,x],na.rm=TRUE)/
    ((2*sum(!is.na(indivlistgeno[indivlistgeno$Year==1990&indivlistgeno$Sex==1,x])))
     +(sum(!is.na(indivlistgeno[indivlistgeno$Year==1990&indivlistgeno$Sex==2,x])))))

#randomly sample from real allele frequencies
simfreq<-sample(datafreq1990,nloci,replace=TRUE)


#sort indivlist
indivlist <- indivlistgeno[order(indivlistgeno$Year),c(1:7)]
colnames(indivlist) <- c( "Year","Indiv", "Category", "Genotyped", "Mom", "Dad", "Sex")

simindivgeno<-indivlist[!duplicated(indivlist$Indiv),]


#separate into moms vs dads vs nestlings
simindivgenoMoms<-
  simindivgeno[(simindivgeno$Category!='nestling' | is.na(simindivgeno$Mom)) & simindivgeno$Sex==2,]
simindivgenoDads<-
  simindivgeno[(simindivgeno$Category!='nestling' | is.na(simindivgeno$Mom)) & simindivgeno$Sex==1,]
simindivgenoNestlings<-
  simindivgeno[simindivgeno$Category=='nestling' & !is.na(simindivgeno$Mom),]

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

####simulate genotypes for nestlings via Mendelian transmission of alleles from parents####

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
  moms.of.sons.geno<-simindivgenoAll[these.moms.of.sons,,drop = FALSE]
  dads.of.daughters.geno<-simindivgenoAll[these.dads.of.daughters,,drop = FALSE]
  dads.of.sons.geno<-simindivgenoAll[these.dads.of.sons,,drop = FALSE]
  
  #check that all parents are genotyped (NOT NA)
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
  simindivgenoAll[these.female.nestlings,8:(nloci+7)] <- female.nestling.geno
  simindivgenoAll[these.male.nestlings,8:(nloci+7)] <- male.nestling.geno
  
  
  #check that the nestlings now have genotypes (not NA)
  stopifnot(all(!is.na(simindivgenoAll[these.female.nestlings,8:(nloci+7)])))
  stopifnot(all(!is.na(simindivgenoAll[these.male.nestlings,8:(nloci+7)])))
  
}

#now we want to have each indiv appear multiple times again
simdataTrue<-merge(indivlist,simindivgenoAll[,c(2,8:(nloci+7))],
                   by.x='Indiv',by.y='Indiv',all.x=TRUE)	


#save
save(simdataTrue,file='simdataTrueZ9May.rdata')

#simdataTrue[c(1:10),c(1:10)]
#colnames(simdataTrue)[7] <- "Sex"

####calculate error in freq estimation due to sampling####

#calculate sample allele freq
#mimic sampling of genotyped indiv by selecting only indivs who actually were genotyped
simdataSample<-simdataTrue[simdataTrue$Genotyped=='Y',]

#get unique indivs in simulated data (all & genotyped)
simdataTrueUnique<-simdataTrue[!duplicated(simdataTrue$Indiv),]
simdataSampleUnique<-simdataSample[!duplicated(simdataSample$Indiv),]
#create data frame to hold simulated allele freqs
simAlleleFreq<-data.frame(Year=integer(),Category=character(),stringsAsFactors=FALSE)

#calculate population (p) and sample allele freq (x)
#and the error in allele freq estimation due to sampling: err = x-p


#parallelize snps
#cores=detectCores() #uncomment these two lines if you want to use more than 4 cores
#cl <- makeCluster(cores[1]-1) #not to overload your computer

#cl <- makeCluster(4) #use 4 cores
#registerDoParallel(cl)

year<-1998
sim<-foreach(i=names(simdataTrue)[8:(nloci+7)],.combine=cbind) %do% {
  #create data frame to hold allele freqs & error
  tmp<-data.frame(Year=rep(year,each=3),Category=c('pt','xt','errT'),
                  stringsAsFactors=FALSE)
  
  frqYr1<-tmp$Year
  frqCat1<-tmp$Category
  
  tmp[frqYr1==year & frqCat1=='pt',3]<-
    sum(simdataTrue[simdataTrue$Year==year,i],na.rm=TRUE)/
    ((2*sum(!is.na(simdataTrue[simdataTrue$Year==year&simdataTrue$Sex==1,i])))
     +(sum(!is.na(simdataTrue[simdataTrue$Year==year&simdataTrue$Sex==2,i]))))
  
  tmp[frqYr1==year & frqCat1=='xt',3]<-
    sum(simdataSample[simdataSample$Year==year,i],na.rm=TRUE)/
    ((2*sum(!is.na(simdataSample[simdataSample$Year==year&simdataSample$Sex==1,i])))
     +(sum(!is.na(simdataSample[simdataSample$Year==year&simdataSample$Sex==2,i]))))
  
  tmp[,3]
}
#stopCluster(cl)

#save the data from this year
save(sim,file=paste("SimAlleleFreqZYr9May_",year,".rdata",sep=''))


#Names for the values we just calculated (year and category/parameter)
simName<-data.frame(Year=rep(year,each=3),Category=c('pt','xt','errT'),
                    stringsAsFactors=FALSE)

#Add the simulation results to the data frame
#err rows are NA b/c we have not calculated the error yet
sim1<-cbind(simName,sim)	
#Add the simulation data to simAlleleFreq (we'll collect the data from all years here)
simAlleleFreq<-rbind(simAlleleFreq,sim1)

#year = 1999
for(year in c(1999:2013)){
  #get moms of sons, dads of sons, and dads of daughters for this year
  
  #get moms of male nestlings born this year
  moms_of_sons<-simdataTrue[simdataTrue$Year==year & simdataTrue$Category=='nestling' & simdataTrue$Sex==1,'Mom']
  #convert list of moms of sons to a data frame
  moms_of_sons<-data.frame(Indiv=moms_of_sons[!is.na(moms_of_sons)],stringsAsFactors=FALSE)
  #collect simulated moms of sons genotypes (including those simulated for ungenotyped indivs) from simdataTrueUnique
  moms_of_sons_geno<-merge(moms_of_sons,simdataTrueUnique[,c(1,8:(nloci+7))],by.x='Indiv',by.y='Indiv',
                           all.x=TRUE)
  
  #simdataTrueUnique[,c(1:8)]
  
  
  #collect simulated moms of sons genotypes (sampled based on real genotyping status) from simdataSampleUnique
  moms_of_sons_genoSample<-merge(moms_of_sons,simdataSampleUnique[,c(1,8:(nloci+7))],by.x='Indiv',
                                 by.y='Indiv',all.x=TRUE)
  #many of the rows in this table are NAs because of the sampling
  
  #don't need moms of daughters as moms do not contribute to daughters
  
  #get dads of male nestlings born this year
  dads_of_sons<-simdataTrue[simdataTrue$Year==year & simdataTrue$Category=='nestling' & simdataTrue$Sex==1,'Dad']
  #convert list of dads of sons to a data frame
  dads_of_sons<-data.frame(Indiv=dads_of_sons[!is.na(dads_of_sons)],stringsAsFactors=FALSE)
  #collect simulated dads of sons genotypes (including those simulated for ungenotyped indivs) from simdataTrueUnique
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
  #collect simulated dads of daughters genotypes (including those simulated for ungenotyped indivs) from simdataTrueUnique
  dads_of_daughters_geno<-merge(dads_of_daughters,simdataTrueUnique[,c(1,8:(nloci+7))],by.x='Indiv',by.y='Indiv',
                                all.x=TRUE)
  #collect simulated dads of daughters genotypes (sampled based on real genotyping status) from simdataSampleUnique
  dads_of_daughters_genoSample<-merge(dads_of_daughters,simdataSampleUnique[,c(1,8:(nloci+7))],by.x='Indiv',
                                      by.y='Indiv',all.x=TRUE)
  #many of the rows in this table are NAs because of the sampling
  
  #parallelize snps
  # cl <- makeCluster(4) #use 4 cores
  # registerDoParallel(cl)
  
  #for each snp
  #snp="1"
  sim<-foreach(i=names(simdataTrue)[8:(nloci+7)],.combine=cbind) %do% {
    #make a data frame to put all these parameters in for each year
    tmp<-data.frame(Year=rep(year,each=69),category=c(
      'pt','xt','errT', 'pt1-pt','xt1-xt', 'errt1-errt',
      
      'pMs','xMs','errMS', 'pMs-pt', 'xMs-xt', 'errMS-errT',
      'pMi','xMi','errMI', 'pMi-pt','xMi-xt','errMI-errT',
      'pMb','xMb','errMB','pMb-pt','xMb-xt','errMB-errT',
      'pMdad','xMdad','errMdad', 
      'pMmom','xMmom','errMmom',
      'pMfam','xMfam','pMfam-pt', 'xMfam-xt', 'errMFAM', 'errMFAM-errT',
      'pMmend','xMmend','errMMEND', 
      
      'pFs','xFs','errFS', 'pFs-pt','xFs-xt','errFS-errT',
      'pFi','xFi','errFI', 'pFi-pt','xFi-xt','errFI-errT',
      'pFb','xFb','errFB', 'pFb-pt','xFb-xt','errFB-errT',
      'pFdad','xFdad','errFdad',
      #take out 'pFmom','xFmom','errFmom',
      'pFfam','xFfam','pFfam-pt','xFfam-xt', 'errFFAM',  'errFFAM-errT',
      'pFmend','xFmend','errFMEND'
    ),stringsAsFactors=FALSE)
    #previously pm = allele frequency of all mothers, pf = allele frequency of all fathers,
    #and similar for xm, xf, errM, errF
    #Now instead we are using pMdad = allele frequency of all fathers of sons,
    #pMmom = allele frequency of all mothers of sons,
    #pFdad = allele frequency of all fathers of daughters,
    #and likewise xMdad, xMmom, errMdad, errMmom, etc.
    
    frqYr1<-tmp$Year
    frqCat1<-tmp$category
    
    #pt is just the mean of all of the simulated data (no sampling)
    tmp[frqYr1==year & frqCat1=='pt',3]<-
      sum(simdataTrue[simdataTrue$Year==year,i],na.rm=TRUE)/
      ((2*sum(!is.na(simdataTrue[simdataTrue$Year==year&simdataTrue$Sex==1,i])))
       +(sum(!is.na(simdataTrue[simdataTrue$Year==year&simdataTrue$Sex==2,i]))))
    
    #xt is the mean of the sampled simulated data
    tmp[frqYr1==year & frqCat1=='xt',3]<-
      sum(simdataSample[simdataSample$Year==year,i],na.rm=TRUE)/
      ((2*sum(!is.na(simdataSample[simdataSample$Year==year&simdataSample$Sex==1,i])))
       +(sum(!is.na(simdataSample[simdataSample$Year==year&simdataSample$Sex==2,i]))))
    
    #ps
    tmp[frqYr1==year & frqCat1=='pMs',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                             simdataTrue$Category=='survivor' & simdataTrue$Sex==1,i])/2
    tmp[frqYr1==year & frqCat1=='pFs',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                             simdataTrue$Category=='survivor' & simdataTrue$Sex==2,i])
    
    #xs
    tmp[frqYr1==year & frqCat1=='xMs',3]<-mean(simdataSample[simdataSample$Year==year & 
                                                               simdataSample$Category=='survivor' & simdataSample$Sex==1,i])/2
    tmp[frqYr1==year & frqCat1=='xFs',3]<-mean(simdataSample[simdataSample$Year==year & 
                                                               simdataSample$Category=='survivor' & simdataSample$Sex==2,i])
    
    #pi
    tmp[frqYr1==year & frqCat1=='pMi',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                             simdataTrue$Category=='immigrant' & simdataTrue$Sex==1,i])/2
    tmp[frqYr1==year & frqCat1=='pFi',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                             simdataTrue$Category=='immigrant' & simdataTrue$Sex==2,i])
    
    #xi
    tmp[frqYr1==year & frqCat1=='xMi',3]<-
      ifelse(is.na(mean(simdataSample[simdataSample$Year==year & 
                                        simdataSample$Category=='immigrant' & simdataSample$Sex==1,i])),0,
             mean(simdataSample[simdataSample$Year==year & 
                                  simdataSample$Category=='immigrant' & simdataSample$Sex==1,i])/2)
    tmp[frqYr1==year & frqCat1=='xFi',3]<-
      ifelse(is.na(mean(simdataSample[simdataSample$Year==year & 
                                        simdataSample$Category=='immigrant' & simdataSample$Sex==2,i])),0,
             mean(simdataSample[simdataSample$Year==year & 
                                  simdataSample$Category=='immigrant' & simdataSample$Sex==2,i]))
    #ifelse here to catch years with no genotyped imms
    #imms are the only category that sometimes is 0 - we checked
    
    #pb
    tmp[frqYr1==year & frqCat1=='pMb',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                             simdataTrue$Category=='nestling' & simdataTrue$Sex==1,i])/2
    tmp[frqYr1==year & frqCat1=='pFb',3]<-mean(simdataTrue[simdataTrue$Year==year & 
                                                             simdataTrue$Category=='nestling' & simdataTrue$Sex==2,i])
    
    #xb
    tmp[frqYr1==year & frqCat1=='xMb',3]<-mean(simdataSample[simdataSample$Year==year &
                                                               simdataSample$Category=='nestling' & simdataSample$Sex==1,i])/2
    tmp[frqYr1==year & frqCat1=='xFb',3]<-mean(simdataSample[simdataSample$Year==year &
                                                               simdataSample$Category=='nestling' & simdataSample$Sex==2,i])
    
    #pMmom & xMmom
    tmp[frqYr1==year & frqCat1=='pMmom',3]<-mean(moms_of_sons_geno[,i],na.rm=TRUE)
    tmp[frqYr1==year & frqCat1=='xMmom',3]<-mean(moms_of_sons_genoSample[,i],na.rm=TRUE)
    
    #pMdad & xMdad
    tmp[frqYr1==year & frqCat1=='pMdad',3]<-mean(dads_of_sons_geno[,i],na.rm=TRUE)/2
    tmp[frqYr1==year & frqCat1=='xMdad',3]<-mean(dads_of_sons_genoSample[,i],na.rm=TRUE)/2
    
    #pFdad & xFdad
    tmp[frqYr1==year & frqCat1=='pFdad',3]<-mean(dads_of_daughters_geno[,i],na.rm=TRUE)/2
    tmp[frqYr1==year & frqCat1=='xFdad',3]<-mean(dads_of_daughters_genoSample[,i],na.rm=TRUE)/2
    
    tmp[,3]
  }
  # stopCluster(cl)
  
  #save p and x results from this year (in case run gets interrupted)
  save(sim,file=paste("SimAlleleFreqZYr9May_",year,".rdata",sep=''))
  
  #categories (parameter names) and years to combine with the results of our calculations
  simName<-data.frame(Year=rep(year,each=69),Category=c(
    'pt','xt','errT', 'pt1-pt', 'xt1-xt', 'errt1-errt',
    
    'pMs','xMs','errMS', 'pMs-pt', 'xMs-xt', 'errMS-errT',
    'pMi','xMi','errMI', 'pMi-pt','xMi-xt','errMI-errT',
    'pMb','xMb','errMB','pMb-pt','xMb-xt','errMB-errT',
    'pMdad','xMdad','errMdad', 
    'pMmom','xMmom','errMmom',
    'pMfam','xMfam','pMfam-pt', 'xMfam-xt', 'errMFAM', 'errMFAM-errT',
    'pMmend','xMmend','errMMEND', 
    
    'pFs','xFs','errFS', 'pFs-pt','xFs-xt','errFS-errT',
    'pFi','xFi','errFI', 'pFi-pt','xFi-xt','errFI-errT',
    'pFb','xFb','errFB', 'pFb-pt','xFb-xt','errFB-errT',
    'pFdad','xFdad','errFdad',
    #take out 'pFmom','xFmom','errFmom',
    'pFfam','xFfam','pFfam-pt','xFfam-xt', 'errFFAM',  'errFFAM-errT',
    'pFmend','xFmend','errFMEND'
  ),stringsAsFactors=FALSE)
  #add the names to the calculation results
  sim1<-cbind(simName,sim)	
  
  #combine with simAlleleFreq which is where we are collecting the results from all the snps
  simAlleleFreq<-rbind(simAlleleFreq,sim1)
}



####calculate error and allele freq differences between each category and the year before####
#err = true error (no hypergeometric error)
frqYr<-simAlleleFreq$Year
frqCat<-simAlleleFreq$Category

#err = true error
simAlleleFreq[frqYr==1998 & frqCat=='errT',c(3:(nloci+2))]<-
  simAlleleFreq[frqYr==1998 & frqCat=='xt',c(3:(nloci+2))]-
  simAlleleFreq[frqYr==1998 & frqCat=='pt',c(3:(nloci+2))]

save(simAlleleFreq,file="simAlleleFreqZ_SR_intermediate9May.rdata")


#calcuate error for each year based on difference b/w p and x
for(year in c(1999:2013)){
  #true error
  simAlleleFreq[frqYr==year & frqCat=='errT',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xt',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pt',c(3:(nloci+2))]	
  
  simAlleleFreq[frqYr==year & frqCat=='errMS',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMs',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pMs',c(3:(nloci+2))]
  simAlleleFreq[frqYr==year & frqCat=='errFS',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFs',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pFs',c(3:(nloci+2))]
  
  simAlleleFreq[frqYr==year & frqCat=='errMI',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMi',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pMi',c(3:(nloci+2))]	
  simAlleleFreq[frqYr==year & frqCat=='errFI',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFi',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pFi',c(3:(nloci+2))]	
  
  simAlleleFreq[frqYr==year & frqCat=='errMB',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMb',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pMb',c(3:(nloci+2))]
  simAlleleFreq[frqYr==year & frqCat=='errFB',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFb',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pFb',c(3:(nloci+2))]
  
  #mom
  simAlleleFreq[frqYr==year & frqCat=='errMmom',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMmom',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pMmom',c(3:(nloci+2))]
  
  #dad
  simAlleleFreq[frqYr==year & frqCat=='errMdad',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMdad',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pMdad',c(3:(nloci+2))]	
  simAlleleFreq[frqYr==year & frqCat=='errFdad',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFdad',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='pFdad',c(3:(nloci+2))]	
  
  
  #allele freq differences
  #true (population) freq (p)
  #this is actually doing pt - pt-1 even thought the category name is pt1-pt
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
  simAlleleFreq[frqYr==year & frqCat=='xt1-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xt',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  
  simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMs',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFs',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  
  simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMi',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFi',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  
  simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMb',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFb',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  
  #true errors
  simAlleleFreq[frqYr==year & frqCat=='errt1-errt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='errT',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
  
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
  
  simAlleleFreq[frqYr==year & frqCat=='xMfam',c(3:(nloci+2))]<-
    0.5*(simAlleleFreq[frqYr==year & frqCat=='xMmom',c(3:(nloci+2))]+
           simAlleleFreq[frqYr==year & frqCat=='xMdad',c(3:(nloci+2))])
  #here we also multiply mom by 0.5 since her contribution makes up only 1/2 of her son's genotype
  
  simAlleleFreq[frqYr==year & frqCat=='xFfam',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFdad',c(3:(nloci+2))]
  #here we don't multiply dad by 0.5 since his contribution makes up 100% of his daughter's genotype
  
  
  simAlleleFreq[frqYr==year & frqCat=='xMmend',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMb',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='xMfam',c(3:(nloci+2))]
  simAlleleFreq[frqYr==year & frqCat=='xFmend',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFb',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==year & frqCat=='xFfam',c(3:(nloci+2))]
  
  simAlleleFreq[frqYr==year & frqCat=='xMfam-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xMfam',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  simAlleleFreq[frqYr==year & frqCat=='xFfam-xt',c(3:(nloci+2))]<-
    simAlleleFreq[frqYr==year & frqCat=='xFfam',c(3:(nloci+2))]-
    simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
  
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

save(simAlleleFreq,file="simAlleleFreqZ9May.rdata")

####calculate variances and covariances####
simVar<-data.frame(Year=rep(c(1999:2013),each=115),Category=rep(c(
  'pt1-pt', 'xt1-xt', 'errt1-errt', 'pt1pterrt1errT',
  'pMs-pt','xMs-xt','errMS-errT','pMspterrMSerrT',
  'pFs-pt','xFs-xt','errFS-errT','pFspterrFSerrT',
  'pMi-pt','xMi-xt','errMI-errT','pMipterrMIerrT',
  'pFi-pt','xFi-xt','errFI-errT','pFipterrFIerrT',
  'pMb-pt','xMb-xt','errMB-errT','pMbpterrMBerrT',
  'pFb-pt','xFb-xt','errFB-errT','pFbpterrFBerrT',
  'pMspMi','xMsxMi','xMserrMI','errMSxMi','errMSerrMI',
  'pMspMb','xMsxMb','xMserrMB','errMSxMb','errMSerrMB',
  'pMipMb','xMixMb','xMierrMB','errMIxMb','errMIerrMB',
  'pFspFi','xFsxFi','xFserrFI','errFSxFi','errFSerrFI',
  'pFspFb','xFsxFb','xFserrFB','errFSxFb','errFSerrFB',
  'pFipFb','xFixFb','xFierrFB','errFIxFb','errFIerrFB',
  'pMspFs','xMsxFs','xMserrFS','errMSxFs','errMSerrFS',
  'pMipFi','xMixFi','xMierrFI','errMIxFi','errMIerrFI',
  'pMbpFb','xMbxFb','xMberrFB','errMBxFb','errMBerrFB',
  'pMspFi','xMsxFi','xMserrFI','errMSxFi','errMSerrFI',
  'pMspFb','xMsxFb','xMserrFB','errMSxFb','errMSerrFB',
  'pMipFb','xMixFb','xMierrFB','errMIxFb','errMIerrFB',
  'pFspMi','xFsxMi','xFserrMI','errFSxMi','errFSerrMI',
  'pFspMb','xFsxMb','xFserrMB','errFSxMb','errFSerrMB',
  'pFipMb','xFixMb','xFierrMB','errFIxMb','errFIerrMB',
  'pMmend','xMmend','errMMEND','pMfam-pt','xMfam-xt','errMFAM-errT',
  'pFmend','xFmend','errFMEND','pFfam-pt','xFfam-xt','errFFAM-errT')
  ,15),stringsAsFactors=FALSE)

bsYr<-simVar$Year
bsCat<-simVar$Category
for(year in c(1999:2013)){
  #total variance, for males and females
  simVar[bsYr==year & bsCat=='pt1-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pt1-pt',c(3:(nloci+2))])^2)
  
  simVar[bsYr==year & bsCat=='xt1-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xt1-xt',c(3:(nloci+2))])^2)
  
  simVar[bsYr==year & bsCat=='errt1-errt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errt1-errt',c(3:(nloci+2))])^2)
  
  simVar[bsYr==year & bsCat=='pt1pterrt1errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pt1-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errt1-errt',c(3:(nloci+2))]))
  
  #variance for each category, for males and females
  #survivors
  simVar[bsYr==year & bsCat=='pMs-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xMs-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errMS-errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='pMspterrMSerrT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pFs-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xFs-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errFS-errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='pFspterrFSerrT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))]))
  
  #imms
  simVar[bsYr==year & bsCat=='pMi-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xMi-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errMI-errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='pMipterrMIerrT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pFi-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xFi-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',3:(nloci+2)])^2)
  simVar[bsYr==year & bsCat=='errFI-errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='pFipterrFIerrT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
  
  #births
  simVar[bsYr==year & bsCat=='pMb-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xMb-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errMB-errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='pMbpterrMBerrT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pFb-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xFb-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errFB-errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='pFbpterrFBerrT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  
  
  #covariance between categories - males
  simVar[bsYr==year & bsCat=='pMspMi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMsxMi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMserrMI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSxMi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSerrMI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pMspMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMsxMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMserrMB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSxMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSerrMB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pMipMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMixMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMierrMB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMIxMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMIerrMB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  
  #covariance between categories - females
  simVar[bsYr==year & bsCat=='pFspFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFsxFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFserrFI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFSxFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFSerrFI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pFspFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFsxFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFserrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFSxFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFSerrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pFipFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFixFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFierrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFIxFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFIerrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  
  
  #covariance between males and females within a category
  #survivors
  simVar[bsYr==year & bsCat=='pMspFs',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMsxFs',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMserrFS',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSxFs',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSerrFS',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))]))
  
  #imms
  simVar[bsYr==year & bsCat=='pMipFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMixFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMierrFI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMIxFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMIerrFI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
  
  #births
  simVar[bsYr==year & bsCat=='pMbpFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMbxFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMberrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMBxFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMBerrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  
  
  #covariance between males and females in different categories
  simVar[bsYr==year & bsCat=='pMspFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMsxFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMserrFI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSxFi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSerrFI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pMspFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMsxFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMserrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSxFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMSerrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pMipFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMixFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xMierrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMIxFb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errMIerrFB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFB-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pFspMi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMi-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFsxMi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFserrMI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFSxMi',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMi-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFSerrMI',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMI-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pFspMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFs-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFsxMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFserrMB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFs-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFSxMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFSerrMB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFS-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  
  simVar[bsYr==year & bsCat=='pFipMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFi-pt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMb-pt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFixMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='xFierrMB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFi-xt',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFIxMb',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMb-xt',c(3:(nloci+2))]))
  simVar[bsYr==year & bsCat=='errFIerrMB',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFI-errT',c(3:(nloci+2))])*
           as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMB-errT',c(3:(nloci+2))]))
  
  
  #variance for mendelian segregation and family size 
  #males
  simVar[bsYr==year & bsCat=='pMmend',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMmend',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xMmend',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMmend',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errMMEND',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMMEND',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='pMfam-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pMfam-pt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xMfam-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xMfam-xt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errMFAM-errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMFAM-errT',c(3:(nloci+2))])^2)
  
  #females
  simVar[bsYr==year & bsCat=='pFmend',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFmend',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xFmend',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFmend',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errFMEND',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFMEND',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='pFfam-pt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pFfam-pt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='xFfam-xt',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xFfam-xt',c(3:(nloci+2))])^2)
  simVar[bsYr==year & bsCat=='errFFAM-errT',3]<-
    mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFFAM-errT',c(3:(nloci+2))])^2)
}

#save output
save(simVar,file="simVarZ9May.rdata")
