## Plotting results of A & Z allele frequency model + simulations to estimate error
## Rose Driscoll and Felix Beaudry
##July 7 2021

setwd('~/Documents/GitHub/ZDropping/alleleFreqModel')

library(plyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scales)
library(dplyr)
library(directlabels)
library(cowplot)
'%ni%' <- function(x,y)!('%in%'(x,y))
library(PNWColors)

#plottheme for manuscript
plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3), 
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"), 
                    plot.background = element_rect(fill = "white"),  
                    axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), 
                    axis.title = element_text(size=7), plot.title = element_text(size=8), 
                    legend.position="right", legend.text = element_text(size=7),
                    legend.title = element_text(size=8), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(0.25,"cm"))


####start calculations for autosomal loci####


load(file='working_files/intermediate_files/indivlistgeno_A.rdata')
indivlist <- indivlistgeno_A[,c(1:7)]

load("working_files/intermediate_files/sampleVar_A.rdata") #samplevar
load("working_files/intermediate_files/simVarA.rdata") #simvar
load("working_files/intermediate_files/allVar_boot_A_w3.4mb.rdata") #bootstrap

col_sample <- sample(length(allVar) -2, 1000, replace=T) + 2

allVar_s <-  allVar[,col_sample]
allVar_q <- allVar[,c(1:2)]

allVar_q$q5 <- apply(allVar_s, 1, function(x) quantile(x,.05,na.rm = T))
allVar_q$q95 <- apply(allVar_s, 1, function(x) quantile(x,.95,na.rm = T))
allVar_q$se <- apply(X=allVar_s,1,function(x) sd(x)/sqrt(length(x)))


## Number of all & Genotyped indivs; proportion of indivs in each category
counts<-ddply(indivlist,.(Year,Category,Sex),summarize,Genotyped=sum(Genotyped=='Y'),
              total=length(Category))
counts['Genotyped']<-2*counts$Genotyped
counts['total']<-2*counts$total

countsAll<-ddply(indivlist,.(Year,Sex),summarize,Genotyped=sum(Genotyped=='Y'),
                 total=length(Category))
countsAll['Genotyped']<-2*countsAll$Genotyped
countsAll['total']<-2*countsAll$total

#proportion of indivs that are in each category and also Genotyped 
propMS<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='survivor' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propFS<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='survivor' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propMI<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='immigrant' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propFI<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='immigrant' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propMB<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='nestling' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propFB<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='nestling' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))

#make tables
prop<-data.frame(propMS=propMS,propFS=propFS,propMI=propMI,propFI=propFI,propMB=propMB,propFB=propFB)
rownames(prop)<-c(2000:2013)

alleleFreqVarAvg<-data.frame(Year=rep(c(2000:2013),each=21),
                             Category=rep(c('MSsq','FSsq','MIsq','FIsq','MBsq','FBsq',
                                            'MSMI','MSMB','MIMB','FSFI','FSFB','FIFB',
                                            'MSFS','MBFB','MIFI',
                                            'MSFI','MSFB','MIFB',
                                            'FSMI','FSMB','FIMB'),14),
                             stringsAsFactors=FALSE)

#calculate final terms for each year, according to formula: (roughly) fraction of population * (change in frequency - error)
#also calculate bootstrap quantiles 
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMS'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMs-xt','avg'] 
                                 - simVar[simVar$Year==x & simVar$Category=='errMS-errT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$Category=='pMspterrMSerrT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMS'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMs-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMS-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMspterrMSerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMS'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMs-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMS-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMspterrMSerrT','q95'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFS'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFs-xt','avg'] 
                                 - simVar[simVar$Year==x & simVar$Category=='errFS-errT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$Category=='pFspterrFSerrT',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFS'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFs-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFS-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFspterrFSerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFS'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFs-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFS-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFspterrFSerrT','q95'])
})



alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMI'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMi-xt','avg'] 
                                 - simVar[simVar$Year==x & simVar$Category=='errMI-errT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$Category=='pMipterrMIerrT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMI'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMi-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMI-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMipterrMIerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMI'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMi-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMI-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMipterrMIerrT','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFI'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFi-xt','avg'] 
                                 - simVar[simVar$Year==x & simVar$Category=='errFI-errT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$Category=='pFipterrFIerrT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFI'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFi-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFI-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFipterrFIerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFI'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFi-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFI-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFipterrFIerrT','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMb-xt','avg'] 
                                 - simVar[simVar$Year==x & simVar$Category=='errMB-errT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$Category=='pMbpterrMBerrT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMb-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMB-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMbpterrMBerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMb-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMB-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMbpterrMBerrT','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='FBsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFb-xt','avg'] 
                                 - simVar[simVar$Year==x & simVar$Category=='errFB-errT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$Category=='pFbpterrFBerrT',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FBsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFb-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFB-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFbpterrFBerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FBsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFb-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFB-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFbpterrFBerrT','q95'])
})




#MM
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxMi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrMI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxMi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrMI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxMi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrMI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxMi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrMI','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxMi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrMI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxMi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrMI','q95'])
})




alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxMb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrMB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxMb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrMB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxMb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrMB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxMb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxMb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrMB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxMb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrMB','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixMb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMierrMB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMIxMb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMIerrMB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIMB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixMb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrMB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxMb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIMB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixMb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrMB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxMb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrMB','q95'])
})





#FF
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxFi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFserrFI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFSxFi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFSerrFI',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxFi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrFI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxFi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxFi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrFI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxFi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrFI','q95'])
})





alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFserrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFSxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFSerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrFB','q95'])
})



alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFixFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFierrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFIxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFIerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFixFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFierrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFixFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFierrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIerrFB','q95'])
})



#MF - same Category
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFS','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFS'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFs','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrFS',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxFs',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrFS',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFS','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFS'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFs','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFS','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFs','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFS','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFS','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFS'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFs','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFS','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFs','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFS','q95'])
})




alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixFi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMierrFI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMIxFi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMIerrFI',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixFi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrFI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxFi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrFI','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixFi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrFI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxFi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrFI','q95'])
})



alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMB']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMbxFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMberrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMBxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMBerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMB']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMbxFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMberrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMBxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMBerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMB']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMbxFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMberrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMBxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMBerrFB','q95'])
})




#MF - different categories
#M then F
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrFI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxFi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrFI',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFI','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFI','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFB','q95'])
})





alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMierrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMIxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMIerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrFB','q95'])
})


#F then M
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxMi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFserrMI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFSxMi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFSerrMI',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxMi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrMI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxMi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrMI','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxMi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrMI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxMi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrMI','q95'])
})



alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxMb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFserrMB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFSxMb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFSerrMB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxMb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrMB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxMb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxMb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrMB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxMb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrMB','q95'])
})




alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFixMb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFierrMB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFIxMb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFIerrMB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIMB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFixMb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFierrMB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIxMb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIMB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFixMb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFierrMB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIxMb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIerrMB','q95'])
})



#Mend and Fam
mendMAvg<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMmend',3] 
                           - simVar[simVar$Year==x & simVar$Category=='errMMEND',3])
})

famMAvg<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMfam-xt',3] 
                           - simVar[simVar$Year==x & simVar$Category=='errMFAM-errT',3])
})

# females
mendFAvg<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFmend',3] 
                           - simVar[simVar$Year==x & simVar$Category=='errFMEND',3])
})

famFAvg<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFfam-xt',3] 
                           - simVar[simVar$Year==x & simVar$Category=='errFFAM-errT',3])
})

newBsq_A <-rbind(data.frame(Year=c(2000:2013),Category=rep('MBsq_mend',14),
                            avg=mendMAvg,stringsAsFactors=FALSE),
                 data.frame(Year=c(2000:2013),Category=rep('MBsq_fam',14),
                            avg=famMAvg,stringsAsFactors=FALSE),
                 data.frame(Year=c(2000:2013),Category=rep('FBsq_mend',14),
                            avg=mendFAvg,stringsAsFactors=FALSE),
                 data.frame(Year=c(2000:2013),Category=rep('FBsq_fam',14),
                            avg=famFAvg,stringsAsFactors=FALSE))

####normalize terms####

#remove irrelevant covariance terms
alleleFreqVarAvg2<-alleleFreqVarAvg[
                                     alleleFreqVarAvg$Category != 'MSMI' 
                                     & alleleFreqVarAvg$Category != 'FSFI'
                                     & alleleFreqVarAvg$Category != 'MSFI'
                                     & alleleFreqVarAvg$Category != 'FSMI'
                                     ,]

# Find sum for each year and calculate proportion for each Category
alleleFreqVarAvg2$sum<-laply(alleleFreqVarAvg2$Year,function(x)   sum(alleleFreqVarAvg2[alleleFreqVarAvg2$Year==x,'avg']))

alleleFreqVarAvg2$prop<-alleleFreqVarAvg2$avg/alleleFreqVarAvg2$sum

alleleFreqVarAvg2$q5_prop <- alleleFreqVarAvg2$q5 / alleleFreqVarAvg2$sum
alleleFreqVarAvg2$q95_prop <- alleleFreqVarAvg2$q95 / alleleFreqVarAvg2$sum

alleleFreqVarAvg1_A <- alleleFreqVarAvg2


####start Z####



#indivlist
load(file='working_files/intermediate_files/indivlistgeno_Z.rdata')
indivlist <- indivlistgeno_Z[,c(1:7)]

load("working_files/intermediate_files/sampleVar_Z.rdata") #samplevar
load("working_files/intermediate_files/simVarZ.rdata") #simvar
load("working_files/intermediate_files/allVar_boot_Z_w3.4mb.rdata") #bootstrap
col_sample <- sample(length(allVar) -2, 1000, replace=T) + 2

allVar_s <-  allVar[,col_sample]
allVar_q <- allVar[,c(1:2)]

allVar_q$q5 <- apply(allVar_s, 1, function(x) quantile(x,.05,na.rm = T))
allVar_q$q95 <- apply(allVar_s, 1, function(x) quantile(x,.95,na.rm = T))
allVar_q$se <- apply(X=allVar_s,1,function(x) sd(x)/sqrt(length(x)))


## Number of all & Genotyped indivs; proportion of indivs in each Category
counts<-ddply(indivlist,.(Year,Category,Sex),summarize,Genotyped=sum(Genotyped=='Y'),
              total=length(Category))
counts[counts$Sex==1,'Genotyped']<-2*counts$Genotyped[counts$Sex==1]
counts[counts$Sex==1,'total']<-2*counts$total[counts$Sex==1]

countsAll<-ddply(indivlist,.(Year,Sex),summarize,Genotyped=sum(Genotyped=='Y'),
                 total=length(Category))
countsAll[countsAll$Sex==1,'Genotyped']<-2*countsAll$Genotyped[countsAll$Sex==1]
countsAll[countsAll$Sex==1,'total']<-2*countsAll$total[countsAll$Sex==1]

propMS<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='survivor' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propFS<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='survivor' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propMI<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='immigrant' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propFI<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='immigrant' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propMB<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='nestling' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
propFB<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$Category=='nestling' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))

prop<-data.frame(propMS=propMS,propFS=propFS,propMI=propMI,propFI=propFI,propMB=propMB,propFB=propFB)
rownames(prop)<-c(2000:2013)

alleleFreqVarAvg<-data.frame(Year=rep(c(2000:2013),each=21),
                             Category=rep(c('MSsq','FSsq','MIsq','FIsq','MBsq','FBsq',
                                            'MSMI','MSMB','MIMB','FSFI','FSFB','FIFB',
                                            'MSFS','MBFB','MIFI',
                                            'MSFI','MSFB','MIFB',
                                            'FSMI','FSMB','FIMB'),14),
                             stringsAsFactors=FALSE)

#square terms
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMS'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMs-xt','avg'] 
                           - simVar[simVar$Year==x & simVar$Category=='errMS-errT',3] 
                           - 2*simVar[simVar$Year==x & simVar$Category=='pMspterrMSerrT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMS'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMs-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMS-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMspterrMSerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMS'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMs-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMS-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMspterrMSerrT','q95'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFS'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFs-xt','avg'] 
                           - simVar[simVar$Year==x & simVar$Category=='errFS-errT',3] 
                           - 2*simVar[simVar$Year==x & simVar$Category=='pFspterrFSerrT',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFS'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFs-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFS-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFspterrFSerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFS'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFs-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFS-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFspterrFSerrT','q95'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMI'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMi-xt','avg'] 
                           - simVar[simVar$Year==x & simVar$Category=='errMI-errT',3] 
                           - 2*simVar[simVar$Year==x & simVar$Category=='pMipterrMIerrT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMI'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMi-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMI-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMipterrMIerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMI'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMi-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMI-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMipterrMIerrT','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFI'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFi-xt','avg'] 
                           - simVar[simVar$Year==x & simVar$Category=='errFI-errT',3] 
                           - 2*simVar[simVar$Year==x & simVar$Category=='pFipterrFIerrT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFI'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFi-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFI-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFipterrFIerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFI'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFi-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFI-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFipterrFIerrT','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMb-xt','avg'] 
                           - simVar[simVar$Year==x & simVar$Category=='errMB-errT',3] 
                           - 2*simVar[simVar$Year==x & simVar$Category=='pMbpterrMBerrT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMb-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMB-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMbpterrMBerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xMb-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMB-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pMbpterrMBerrT','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='FBsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFb-xt','avg'] 
                           - simVar[simVar$Year==x & simVar$Category=='errFB-errT',3] 
                           - 2*simVar[simVar$Year==x & simVar$Category=='pFbpterrFBerrT',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FBsq','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFb-xt','q5'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFB-errT','q5'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFbpterrFBerrT','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FBsq','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(allVar_q[allVar_q$Year==x & allVar_q$Category=='xFb-xt','q95'] 
                           - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFB-errT','q95'] 
                           - 2*allVar_q[allVar_q$Year==x & allVar_q$Category=='pFbpterrFBerrT','q95'])
})




#MM
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxMi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrMI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxMi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrMI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxMi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrMI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxMi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrMI','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxMi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrMI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxMi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrMI','q95'])
})




alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxMb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrMB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxMb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrMB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxMb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrMB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxMb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxMb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrMB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxMb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrMB','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixMb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMierrMB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMIxMb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMIerrMB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIMB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixMb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrMB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxMb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIMB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixMb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrMB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxMb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrMB','q95'])
})





#FF
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxFi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFserrFI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFSxFi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFSerrFI',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxFi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrFI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxFi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxFi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrFI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxFi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrFI','q95'])
})





alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFserrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFSxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFSerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrFB','q95'])
})



alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFixFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFierrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFIxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFIerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFixFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFierrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFixFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFierrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIerrFB','q95'])
})



#MF - same Category
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFS','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFS'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFs','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrFS',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxFs',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrFS',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFS','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFS'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFs','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFS','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFs','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFS','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFS','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFS'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFs','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFS','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFs','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFS','q95'])
})




alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixFi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMierrFI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMIxFi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMIerrFI',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixFi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrFI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxFi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrFI','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixFi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrFI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxFi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrFI','q95'])
})



alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMB']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMbxFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMberrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMBxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMBerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMB']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMbxFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMberrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMBxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMBerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMB']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMbxFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMberrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMBxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMBerrFB','q95'])
})




#MF - different categories
#M then F
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrFI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxFi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrFI',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFI','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFI','q95'])
})






alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMserrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMSxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMSerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMS']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMsxFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMserrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMSerrFB','q95'])
})





alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixFb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xMierrFB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errMIxFb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errMIerrFB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixFb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrFB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxFb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrFB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propMI']*prop[yr,'propFB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xMixFb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxMierrFB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIxFb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errMIerrFB','q95'])
})


#F then M
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxMi','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFserrMI',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFSxMi',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFSerrMI',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMI','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxMi','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrMI','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxMi','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrMI','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMI','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMI'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxMi','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrMI','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxMi','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrMI','q95'])
})



alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxMb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFserrMB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFSxMb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFSerrMB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxMb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrMB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxMb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFS']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFsxMb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFserrMB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSxMb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFSerrMB','q95'])
})




alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFixMb','avg'] 
     - simVar[simVar$Year==x & simVar$Category=='xFierrMB',3] 
     - simVar[simVar$Year==x & simVar$Category=='errFIxMb',3] 
     + simVar[simVar$Year==x & simVar$Category=='errFIerrMB',3])
})
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIMB','q5']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFixMb','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFierrMB','q5'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIxMb','q5'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIerrMB','q5'])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIMB','q95']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(prop[yr,'propFI']*prop[yr,'propMB'])*
    (allVar_q[allVar_q$Year==x & allVar_q$Category=='xFixMb','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='sxFierrMB','q95'] 
     - allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIxMb','q95'] 
     + allVar_q[allVar_q$Year==x & allVar_q$Category=='errFIerrMB','q95'])
})


#Mend and Fam
mendMAvg<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMmend',3] 
                           - simVar[simVar$Year==x & simVar$Category=='errMMEND',3])
})


famMAvg<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propMB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMfam-xt',3] 
                           - simVar[simVar$Year==x & simVar$Category=='errMFAM-errT',3])
})

# females
mendFAvg<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFmend',3] 
                           - simVar[simVar$Year==x & simVar$Category=='errFMEND',3])
})

famFAvg<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  ((prop[yr,'propFB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFfam-xt',3] 
                           - simVar[simVar$Year==x & simVar$Category=='errFFAM-errT',3])
})

newBsq_Z <-rbind(data.frame(Year=c(2000:2013),Category=rep('MBsq_mend',14),
                            avg=mendMAvg,stringsAsFactors=FALSE),
                 data.frame(Year=c(2000:2013),Category=rep('MBsq_fam',14),
                            avg=famMAvg,stringsAsFactors=FALSE),
                 data.frame(Year=c(2000:2013),Category=rep('FBsq_mend',14),
                            avg=mendFAvg,stringsAsFactors=FALSE),
                 data.frame(Year=c(2000:2013),Category=rep('FBsq_fam',14),
                            avg=famFAvg,stringsAsFactors=FALSE))

####normalize Z terms####

alleleFreqVarAvg2<-alleleFreqVarAvg[
                              
                                     alleleFreqVarAvg$Category != 'MSMI' 
                                    & alleleFreqVarAvg$Category != 'FSFI'
                                    & alleleFreqVarAvg$Category != 'MSFI'
                                    & alleleFreqVarAvg$Category != 'FSMI'
                                    
                                    # Remove immigrant covariance and mom-daughter Z inheritance
                                    
                                    & alleleFreqVarAvg$Category != 'FIFB'
                                    & alleleFreqVarAvg$Category != 'FSFB',]


# Find sum for each year and calculate proportion for each Category
alleleFreqVarAvg2$sum<-laply(alleleFreqVarAvg2$Year,function(x)   sum(alleleFreqVarAvg2[alleleFreqVarAvg2$Year==x,'avg']))

alleleFreqVarAvg2$prop<-alleleFreqVarAvg2$avg/alleleFreqVarAvg2$sum
alleleFreqVarAvg2$q5_prop <- alleleFreqVarAvg2$q5 / alleleFreqVarAvg2$sum
alleleFreqVarAvg2$q95_prop <- alleleFreqVarAvg2$q95 / alleleFreqVarAvg2$sum

alleleFreqVarAvg1_Z <- alleleFreqVarAvg2

####combine A and Z####
#stack
AZ_AFVA<-rbind.data.frame(
  cbind.data.frame(alleleFreqVarAvg1_A,"chrom"="A"),
  cbind.data.frame(alleleFreqVarAvg1_Z,"chrom"="Z")
)
#merge
alleleFreqVarAvg1_AZ <- left_join(alleleFreqVarAvg1_A,alleleFreqVarAvg1_Z,by=c("Year"="Year","Category"="Category"))

####categories and sexes####
detach(package:plyr)

#set sex label & factor
AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("FSsq","FIsq","FBsq","FIFB","FSFB")] <- "Female"
AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("MSsq","MIsq","MBsq","MIMB","MSMB")] <- "Male"
AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("FIMB","FSMB","MIFB","MSFB","MBFB","MSFS","MIFI")] <- "Cov(F,M)"
AZ_AFVA$sexCat <- factor(AZ_AFVA$SexCat,levels=c("Female" ,"Male"    ,   "Cov(F,M)"))

AZ_AFVA_sex <- 
  AZ_AFVA %>% group_by(chrom,Year) %>% 
  mutate(sex_sum=sum(avg)) %>% 
  group_by(chrom,Year,SexCat) %>% 
  mutate(sexCat_sum=sum(avg)/sex_sum) 

AZ_AFVA_sex_only <- unique(AZ_AFVA_sex[,c(1,10,12:14)])
AZ_AFVA_sex_only$SexCat = factor(AZ_AFVA_sex_only$sexCat, levels=c('Cov(F,M)','Female','Male'))



#set factor order 
AZ_AFVA$faccat <- factor(AZ_AFVA$Category,
                         levels = 
                           c("MSsq", "MIsq","MBsq",
                             "MSMB" , "MIMB"  ,"FBsq",
                             "FSsq"  , "FIsq",
                             "FSFB"   , "FIFB"    , 
                             "MSFS"   ,   "MBFB"  ,     "MIFI"   , 
                             "MSFB"    ,   "MIFB"   ,   "FSMB"   ,   "FIMB"     
                           ))

#set demographic group label & factor
AZ_AFVA$Supercategory<-factor(ifelse(AZ_AFVA$Category=='MIFI' | AZ_AFVA$Category=='MIsq' |  AZ_AFVA$Category=='FIsq','Immigrant',
                                     ifelse(AZ_AFVA$Category=='MIMB' | AZ_AFVA$Category=='FIFB' | AZ_AFVA$Category=='MIFB' |  AZ_AFVA$Category=='FIMB', 'Cov(I,B)',
                                            ifelse( AZ_AFVA$Category=='MBFB' | AZ_AFVA$Category=='MBsq' | AZ_AFVA$Category=='FBsq', 'Birth',
                                                    ifelse(AZ_AFVA$Category=='MSMB' | AZ_AFVA$Category=='FSFB' | AZ_AFVA$Category=='MSFB' | AZ_AFVA$Category=='FSMB', 'Cov(S,B)',
                                                           ifelse(AZ_AFVA$Category=='MSsq' |AZ_AFVA$Category=='FSsq' | AZ_AFVA$Category=='MSFS', 'Survivor', NA))))),
                              levels = c('Immigrant', 'Cov(I,B)', 'Birth','Cov(S,B)','Survivor'))

AZ_AFVA$supercategory <- factor(AZ_AFVA$Supercategory,levels=c("Survivor","Cov(S,B)","Birth","Cov(I,B)","Immigrant"))

AZ_AFVA$Category2 <- AZ_AFVA$Category
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("FSsq","MSsq","MSFS")] <- "S"
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("FIsq","MIsq","MIFI")] <- "I"
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("FBsq","MBsq","MBFB")] <- "B"
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("MIMB","FIFB")] <- "IB"
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("MSMB","FSFB")] <- "SB"

#set sex covariances
AZ_AFVA$Category4 <- NA
AZ_AFVA$Category4[AZ_AFVA$Category %in% c("FIMB","FSMB")] <- "FM"
AZ_AFVA$Category4[AZ_AFVA$Category %in% c("MIFB","MSFB")] <- "MF"


AZ_AFVA_cat <- 
  AZ_AFVA %>% group_by(chrom,Year) %>% 
  mutate(cat_sum=sum(avg)) %>% 
  group_by(chrom,Year,supercategory) %>% 
  mutate(catCat_sum=sum(avg)/cat_sum) 


AZ_AFVA_cat_only <- unique(AZ_AFVA_cat[,c(1,10,15,18:19)])
AZ_AFVA_cat_only$supercategory = factor(AZ_AFVA_cat_only$supercategory, levels=c('Cov(I,B)','Immigrant','Cov(S,B)','Birth','Survivor'))

fills_to_use_cat <-c(pnw_palettes$Bay[1,1],pnw_palettes$Bay[1,2],pnw_palettes$Bay[1,5], pnw_palettes$Bay[1,3],pnw_palettes$Bay[1,4])

####plot fig. 4B####
varp_title <- expression(atop("Proportional Contribution",paste("to ",Delta, "p variance")))


sex_stack_plot <- 
  ggplot() + 
  geom_hline(yintercept = 0,alpha=0.5)+
  geom_hline(yintercept = 1,alpha=0.5)+
  
  geom_bar(data=AZ_AFVA_sex_only, aes(x=Year, y=sexCat_sum,fill=SexCat,group=SexCat), stat='identity',alpha=0.75) +
  
  facet_grid(chrom~.,scales="free",space="free") + 
  labs(y=varp_title,x="Year") +
  theme_bw(base_size = 8) + 
  theme(strip.background =element_rect(fill="white")) +
  scale_fill_manual("Category",values=c("mediumpurple","indianred1","cornflowerblue")) + 
  guides(color=FALSE)  + 
  theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_y_continuous(breaks=c(0,0.5,1))

min(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Male" & AZ_AFVA_sex_only$chrom == "A"])
max(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Male" & AZ_AFVA_sex_only$chrom == "A"])

min(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Female" & AZ_AFVA_sex_only$chrom == "A"])
max(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Female" & AZ_AFVA_sex_only$chrom == "A"])

min(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Male" & AZ_AFVA_sex_only$chrom == "Z"])
max(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Male" & AZ_AFVA_sex_only$chrom == "Z"])

min(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Female" & AZ_AFVA_sex_only$chrom == "Z"])
max(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Female" & AZ_AFVA_sex_only$chrom == "Z"])


cat_stack_plot <- 
  ggplot() + 
  geom_hline(yintercept = 0,alpha=0.5)+
  
  geom_hline(yintercept = 1,alpha=0.5)+
  
  geom_bar(data=AZ_AFVA_cat_only, aes(x=Year, y=catCat_sum,fill=supercategory,group=supercategory), stat='identity',alpha=0.75) +
  
  facet_grid(chrom~.,scales="free",space="free") + 
  labs(y=varp_title,x="Year") +
  theme_bw(base_size = 8) + 
  theme(strip.background =element_rect(fill="white")) +
  scale_fill_manual("Category",values=fills_to_use_cat) + 
  guides(color=FALSE)  + 
  theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_y_continuous(breaks=c(0,0.5,1))

rep_sum_A <- AZ_AFVA_cat_only[AZ_AFVA_cat_only$supercategory %in% c("Survivor","Birth","Cov(S,B)") & AZ_AFVA_cat_only$chrom == "A",] %>% 
  group_by(chrom,Year) %>% 
  summarize(rep_sum=sum( catCat_sum)) 

rep_sum_Z <- AZ_AFVA_cat_only[AZ_AFVA_cat_only$supercategory %in% c("Survivor","Birth","Cov(S,B)") & AZ_AFVA_cat_only$chrom == "Z",] %>% 
  group_by(chrom,Year) %>% 
  summarize(rep_sum=sum( catCat_sum)) 


#quartz(width=6, height=4,dpi=300)
pdf(paste("fig_4BC.pdf",sep=''),width=4.5,height=3)

plot_grid(sex_stack_plot,  cat_stack_plot,  labels = c('B', 'C'), ncol = 1, align = 'v',axis='rl')
dev.off()




####label supercategories for stack####



#plotting parameters
Avarp_title <- expression(paste("Autosomal Variance ",Delta, "p"))
Zvarp_title <- expression(paste("Z Variance ",Delta, "p"))

fills_to_use <-c(pnw_palettes$Bay[1,3],pnw_palettes$Winter[1,2],"#DD4124", pnw_palettes$Bay[1,2],pnw_palettes$Bay[1,1],pnw_palettes$Winter[1,4],"#E77A66",pnw_palettes$Bay[1,4],pnw_palettes$Bay[1,5])

####Plot Fig. 5####

AZ_AFVA$year_adj <- AZ_AFVA$Year 
AZ_AFVA$year_adj[AZ_AFVA$chrom == "A"] <- AZ_AFVA$year_adj[AZ_AFVA$chrom == "A"] + 0.25

varp_title2 <- expression(paste("Proportional Contribution to ",Delta, "p variance"))

pdf(paste("fig_5.pdf",sep=''),width=5.5,height=4)

ggplot(data=AZ_AFVA, aes(x=year_adj, y=prop)) + 
  geom_hline(yintercept = 0,alpha=0.5,size=0.25)+
  
  geom_errorbar(aes(ymin=q5_prop, ymax=q95_prop,color=Category2,linetype=chrom), width=.1,alpha=0.5,size=0.15) +
#  geom_line(aes(linetype=chrom,color=Category2),alpha=0.5) + 
  
  geom_point(aes(color=Category2,shape=chrom),size=0.75) + 
 guides(color=FALSE) +
  facet_grid(  supercategory ~sexCat,scales="free")+
  
  theme_bw(base_size = 8) + 
  scale_x_continuous(breaks=c(2000,2005,2010),limits = c(1999,2015))+  
  scale_color_manual(values=fills_to_use)+  
  scale_shape_manual(values=c(16,1)) +
 # scale_linetype_manual(values=c("solid","dashed")) +
  labs(y=varp_title2,x="Year",color="Category",linetype="") + 
  theme(strip.background =element_rect(fill="white")) +
  geom_dl(aes(label = Category4,color=Category2), method = list("last.qp",cex = 0.5,dl.trans(x = x + .2))) +
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + 
  coord_cartesian(ylim = c(-0.25, 0.5))
  #ylim(-0.2,0.5)
dev.off()


min(AZ_AFVA$prop[AZ_AFVA$supercategory == "Immigrant"])
max(AZ_AFVA$prop[AZ_AFVA$supercategory == "Immigrant"])

 min(AZ_AFVA$prop[AZ_AFVA$supercategory == "Birth"])
max(AZ_AFVA$prop[AZ_AFVA$supercategory == "Birth"])

min(AZ_AFVA$prop[AZ_AFVA$supercategory == "Survivor"])
max(AZ_AFVA$prop[AZ_AFVA$supercategory == "Survivor"])

####label supercategories for merge####

names(alleleFreqVarAvg1_AZ) <- 
c("Year"   ,    "Category" ,"avg.A"   ,   "q5.A"   ,    "q95.A"    ,  "sum.A"  ,   
"prop.A"  ,   "q5_prop.A",  "q95_prop.A","avg.Z"    ,  "q5.Z"     ,  "q95.Z" ,    
 "sum.Z"      ,"prop.Z"  ,   "q5_prop.Z",  "q95_prop.Z")

#rename in full
alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category == "MSsq"] <- "Male Survivor"
alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MIsq"] <- "Male Immigrant"
alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MBsq"] <- "Male Birth"
alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category == "FSsq"] <- "Female Survivor"
alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category == "FIsq" ] <- "Female Immigrant"
alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category == "FBsq" ] <- "Female Birth"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MSMB" ] <- "Cov(Ms,Mb)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MIMB" ] <- "Cov(Mi,Mb)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "FSFB"] <- "Cov(Fs,Fb)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "FIFB"] <- "Cov(Fi,Fb)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MSFS" ]<- "Cov(Ms,Fs)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MBFB" ]<- "Cov(Mb,Fb)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MIFI"]<- "Cov(Mi,Fi)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MSFB"]<- "Cov(Ms,Fb)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "MIFB"]<- "Cov(Mi,Fb)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "FSMB" ]<- "Cov(Fs,Mb)"
 alleleFreqVarAvg1_AZ$fullCat[alleleFreqVarAvg1_AZ$Category ==  "FIMB"]<- "Cov(Fi,Mb)"

#set factors
 alleleFreqVarAvg1_AZ$FullCat <- factor(alleleFreqVarAvg1_AZ$fullCat, 
                                levels=c(
 "Male Survivor"   , "Female Survivor" , "Cov(Ms,Fs)"   , "Cov(Mb,Fb)"    ,   "Cov(Ms,Fb)"    ,
 "Male Birth"    ,   "Female Birth"   ,   "Cov(Ms,Mb)"     ,  "Cov(Fs,Mb)"    ,  "Cov(Mi,Fb)"   ,
 "Male Immigrant"  , "Female Immigrant",     "Cov(Mi,Fi)"    ,   "Cov(Mi,Mb)"    ,     
        "Cov(Fi,Mb)"  ,"Cov(Fs,Fb)","Cov(Fi,Fb)"))

 #Fig S6
 pdf(paste("fig_S7.pdf",sep=''),width=6,height=4)
 
 ggplot(alleleFreqVarAvg1_AZ[alleleFreqVarAvg1_AZ$FullCat %ni% c("Cov(Fs,Fb)","Cov(Fi,Fb)"),],aes(x=prop.A,y=prop.Z)) +   
  geom_abline(intercept = 0,slope=(4/3),alpha=0.5,linetype="dotted",size=0.25) +
  geom_abline(intercept = 0,slope=1,alpha=0.5,size=0.25) +
  geom_abline(intercept = 0,slope=(2/3),alpha=0.5,linetype="dashed",size=0.25) +
  
  facet_wrap(~ FullCat,scales="free",ncol=5) + 
  geom_smooth(method="lm",color="black",size=0.25) +  
  geom_point(size=0.25) + 
  #geom_text(aes(color=as.factor(Year),label=as.factor(Year))) + 
  labs(x=Avarp_title,y=Zvarp_title,color="Year") + 
  guides(color=FALSE) +
  theme_bw(base_size = 8) + 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
  theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
theme(strip.background =element_rect(fill="white")) 

 dev.off()
 



####Mend and Fam####
##Autosomal
#re-calculate proportions without births
alleleFreqVarAvg1_newBsq_A <- 
rbind.data.frame(
  alleleFreqVarAvg1_A[alleleFreqVarAvg1_A$Category %ni% c("MBsq","FBsq"),c(1:3)],
  newBsq_A
)
alleleFreqVarAvg1_newBsq_A$sum<-lapply(alleleFreqVarAvg1_newBsq_A$Year,function(x)   sum(alleleFreqVarAvg1_newBsq_A[alleleFreqVarAvg1_newBsq_A$Year==x,'avg']))
alleleFreqVarAvg1_newBsq_A$avg <- as.numeric(alleleFreqVarAvg1_newBsq_A$avg)
alleleFreqVarAvg1_newBsq_A$sum <- as.numeric(alleleFreqVarAvg1_newBsq_A$sum)

alleleFreqVarAvg1_newBsq_A$prop<-alleleFreqVarAvg1_newBsq_A$avg/alleleFreqVarAvg1_newBsq_A$sum

#re-add births for comparison
alleleFreqVarAvg1_newBsq_A <- 
  rbind.data.frame(
    alleleFreqVarAvg1_A[alleleFreqVarAvg1_A$Category %in% c("MBsq","FBsq"),c(1,2,3,6,7)],
    alleleFreqVarAvg1_newBsq_A
  )
alleleFreqVarAvg1_newBsq_A$chrom <- "Autosomal"

#Z
alleleFreqVarAvg1_newBsq_Z <- 
  rbind.data.frame(
    alleleFreqVarAvg1_Z[alleleFreqVarAvg1_Z$Category %ni% c("MBsq","FBsq"),c(1:3)],
    newBsq_Z
  )

alleleFreqVarAvg1_newBsq_Z$sum<-lapply(alleleFreqVarAvg1_newBsq_Z$Year,function(x)   sum(alleleFreqVarAvg1_newBsq_Z[alleleFreqVarAvg1_newBsq_Z$Year==x,'avg']))
alleleFreqVarAvg1_newBsq_Z$avg <- as.numeric(alleleFreqVarAvg1_newBsq_Z$avg)
alleleFreqVarAvg1_newBsq_Z$sum <- as.numeric(alleleFreqVarAvg1_newBsq_Z$sum)

alleleFreqVarAvg1_newBsq_Z$prop<-alleleFreqVarAvg1_newBsq_Z$avg/alleleFreqVarAvg1_newBsq_Z$sum

alleleFreqVarAvg1_newBsq_Z <- 
  rbind.data.frame(
    alleleFreqVarAvg1_Z[alleleFreqVarAvg1_Z$Category %in% c("MBsq","FBsq"),c(1,2,3,6,7)],
    alleleFreqVarAvg1_newBsq_Z
  )
alleleFreqVarAvg1_newBsq_Z$chrom <- "Z"

alleleFreqVarAvg1_newBsq <- rbind.data.frame(
  alleleFreqVarAvg1_newBsq_Z,
  alleleFreqVarAvg1_newBsq_A
)

newBsq <- alleleFreqVarAvg1_newBsq[alleleFreqVarAvg1_newBsq$Category %in% c("FBsq_fam","FBsq_mend","MBsq_fam","MBsq_mend","MBsq","FBsq"),]

newBsq$sexCat[newBsq$Category %in% c("FBsq_fam","FBsq_mend","FBsq")] <- "Female"
newBsq$sexCat[newBsq$Category %in% c("MBsq_fam","MBsq_mend","MBsq")] <- "Male"

newBsq$bCat[newBsq$Category %in% c("FBsq_mend","MBsq_mend")] <- "Mendelian Noise"
newBsq$bCat[newBsq$Category %in% c("FBsq_fam","MBsq_fam")] <- "Family Size"
newBsq$bCat[newBsq$Category %in% c("MBsq","FBsq")] <- "Birth Total"

#varp_birth_title <- expression(paste("Variance ",Delta, "p Birth"))

#fig S8
#mendFam_year_plot <- 
pdf(paste("fig_S9.pdf",sep=''),width=6,height=4)

ggplot() + 
  geom_hline(yintercept = 0,alpha=0.2)+
 # geom_line(data=newBsq[newBsq$bCat == "Birth Total",], aes(x=Year, y=prop,linetype=bCat),alpha=0.2,linetype="longdash") +
  geom_line(data=newBsq[newBsq$bCat %in% c("Family Size","Mendelian Noise"),], aes(x=Year, y=prop,color=sexCat,group=bCat),position="stack",alpha=0.2) +
  geom_point(data=newBsq[newBsq$bCat %in% c("Family Size","Mendelian Noise"),], aes(x=Year, y=prop,color=sexCat,shape=bCat),position="stack") +
  facet_grid(chrom~sexCat) + 
  labs(y=varp_title,x="Year",linetype="Birth Term",color="Sex",shape="Birth Term") +
  theme_bw(base_size = 8) + 
  theme(strip.background =element_rect(fill="white")) +
   guides(color=FALSE)  + 
  scale_color_manual(values=c("indianred1","cornflowerblue" )) + 
  theme( panel.grid.minor = element_blank()) +
  scale_shape_manual(values=c(3,4)) +
  theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) 

dev.off()


newBsq_ZA <- left_join(
  newBsq[newBsq$chrom == "Z",c(1:5)],
  newBsq[newBsq$chrom == "Autosomal",c(1:5)], by=c("Year"="Year","Category"="Category")
)
names(newBsq_ZA) <- c("Year","Category","avg.Z","sum.Z","prop.Z","avg.A","sum.A","prop.A")

newBsq_ZA$SexCat[newBsq_ZA$Category == "FBsq_fam"]  <- "Female"
newBsq_ZA$Birth[newBsq_ZA$Category == "FBsq_fam"] <- "Family Size"

newBsq_ZA$SexCat[newBsq_ZA$Category == "FBsq_mend"] <- "Female"
newBsq_ZA$Birth[newBsq_ZA$Category == "FBsq_mend"] <- "Mendelian Noise"

newBsq_ZA$SexCat[newBsq_ZA$Category == "MBsq_fam"]<- "Male"
newBsq_ZA$Birth[newBsq_ZA$Category == "MBsq_fam"] <- "Family Size"

newBsq_ZA$SexCat[newBsq_ZA$Category == "MBsq_mend"] <- "Male"
newBsq_ZA$Birth[newBsq_ZA$Category == "MBsq_mend"] <- "Mendelian Noise"

#mendFam_ZA_plot <- 
newBsq_ZA <- newBsq_ZA[!is.na(newBsq_ZA$Birth),]
#fig s9  
pdf(paste("fig_S10.pdf",sep=''),width=6,height=4)

ggplot(newBsq_ZA,aes(x=prop.A,y=prop.Z)) + 
  geom_abline(intercept = 0,slope=(4/3),alpha=0.5,linetype="dotted") +
  
  geom_abline(intercept = 0,slope=1,alpha=0.5,linetype="dashed") +
  geom_abline(intercept = 0,slope=(2/3),alpha=0.5,linetype="dotted") +
  
  facet_grid(Birth~ SexCat,scales="free") + 
  geom_smooth(method="lm",color="black") +  
  geom_point(aes(color=as.factor(SexCat),shape=Birth)) + 
  #geom_text(aes(color=as.factor(Year),label=as.factor(Year))) + 
  labs(x=Avarp_title,y=Zvarp_title,color="Sex",shape="Birth Term") + 
  guides() +
  scale_shape_manual(values=c(3,4)) + 
  theme_bw(base_size = 8) + 
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
  theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(strip.background =element_rect(fill="white")) +
  scale_color_manual(values=c("indianred1","cornflowerblue" )) +
  theme(panel.spacing = unit(0.3, "cm")) 
 
dev.off()




