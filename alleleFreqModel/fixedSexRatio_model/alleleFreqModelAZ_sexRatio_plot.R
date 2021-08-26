## Plotting results of A and Z allele frequency change partitioning model with fixed (50:50) sex ratio; prints figures equivalent to Fig 4-5 and S6-10 
## last update: Aug 26th 2021
## Rose Driscoll and Felix Beaudry


## Setup

library(plyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scales)
library(dplyr)
library(directlabels)

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

setwd('~/Google Drive/Research/Data2/fsj/ZDropping_all')

#indivlist
load("simindivFIXmin2obs.rdata")
ped<-read.table('FSJpedgeno_Zsexlinked.ped',header=FALSE,sep=' ',stringsAsFactors=FALSE)
pedinfo <- ped[,1:5]
colnames(pedinfo) <- c("Family", "USFWS", "Dad", "Mom", "Sex")
indivlist <- merge(simindivFIXmin2obs[,1:6],pedinfo[,c(2,5)],by='USFWS')
indivlist <- indivlist[order(indivlist$Year),]

#samplevar
load("sampleVarA_06May2021_rose.rdata")

#simvar
load("simVarA_07Apr2021_rose.rdata")


## Number of all & genotyped indivs; proportion of indivs in each category
#need to multiply by 2 for everyone here (not just males)

counts<-ddply(indivlist,.(Year,category,Sex),summarize,genotyped=sum(genotyped=='Y'),
              total=length(category))
counts[counts$Sex==1,'genotyped']<-2*counts$genotyped[counts$Sex==1]
counts[counts$Sex==2,'genotyped']<-2*counts$genotyped[counts$Sex==2]
counts[counts$Sex==1,'total']<-2*counts$total[counts$Sex==1]
counts[counts$Sex==2,'total']<-2*counts$total[counts$Sex==2]
countsAll<-ddply(indivlist,.(Year,Sex),summarize,genotyped=sum(genotyped=='Y'),
                 total=length(category))
countsAll[countsAll$Sex==1,'genotyped']<-2*countsAll$genotyped[countsAll$Sex==1]
countsAll[countsAll$Sex==2,'genotyped']<-2*countsAll$genotyped[countsAll$Sex==2]
countsAll[countsAll$Sex==1,'total']<-2*countsAll$total[countsAll$Sex==1]
countsAll[countsAll$Sex==2,'total']<-2*countsAll$total[countsAll$Sex==2]

propMS<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='survivor' & counts$Sex==1,'total']/countsAll[countsAll$Year==x & countsAll$Sex==1,'total'])
propFS<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='survivor' & counts$Sex==2,'total']/countsAll[countsAll$Year==x & countsAll$Sex==2,'total'])
propMI<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='immigrant' & counts$Sex==1,'total']/countsAll[countsAll$Year==x & countsAll$Sex==1,'total'])
propFI<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='immigrant' & counts$Sex==2,'total']/countsAll[countsAll$Year==x & countsAll$Sex==2,'total'])
propMB<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='nestling' & counts$Sex==1,'total']/countsAll[countsAll$Year==x & countsAll$Sex==1,'total'])
propFB<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='nestling' & counts$Sex==2,'total']/countsAll[countsAll$Year==x & countsAll$Sex==2,'total'])
# I think this is proportion of indivs that are in each category and also genotyped? 
# so it is out of total # indivs for that year
prop<-data.frame(propMS=propMS,propFS=propFS,propMI=propMI,propFI=propFI,propMB=propMB,propFB=propFB)
rownames(prop)<-c(2000:2013)
#rowSums(prop)
#Adds to 2 - is that correct?
#Yes because MS+MI+MB should add to 1 and FS+FI+FB should add to 1


alleleFreqVarAvg<-data.frame(Year=rep(c(2000:2013),each=21),
                             Category=rep(c('MSsq','FSsq','MIsq','FIsq','MBsq','FBsq',
                                            'MSMI','MSMB','MIMB','FSFI','FSFB','FIFB',
                                            'MSFS','MBFB','MIFI',
                                            'MSFI','MSFB','MIFB',
                                            'FSMI','FSMB','FIMB'),14),
                             stringsAsFactors=FALSE)

#square terms
x=2000
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/2)*prop[yr,'propMS'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMs-xMt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errMS-errMT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pMspMterrMSerrMT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/2)*prop[yr,'propFS'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFs-xFt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errFS-errFT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pFspFterrFSerrFT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/2)*prop[yr,'propMI'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMi-xMt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errMI-errMT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pMipMterrMIerrMT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/2)*prop[yr,'propFI'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFi-xFt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errFI-errFT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pFipFterrFIerrFT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/2)*prop[yr,'propMB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMb-xMt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errMB-errMT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pMbpMterrMBerrMT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FBsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/2)*prop[yr,'propFB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFb-xFt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errFB-errFT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pFbpFterrFBerrFT',3])
})

#MM
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/2)^2)*(prop[yr,'propMS']*prop[yr,'propMI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxMi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrMI',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxMi',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrMI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/2)^2)*(prop[yr,'propMS']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxMb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrMB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxMb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrMB',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/2)^2)*(prop[yr,'propMI']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixMb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMierrMB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMIxMb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMIerrMB',3])
})

#FF
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/2)^2)*(prop[yr,'propFS']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxFi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFserrFI',3] 
     - simVar[simVar$Year==x & simVar$category=='errFSxFi',3] 
     + simVar[simVar$Year==x & simVar$category=='errFSerrFI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/2)^2)*(prop[yr,'propFS']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFserrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errFSxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errFSerrFB',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/2)^2)*(prop[yr,'propFI']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFixFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFierrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errFIxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errFIerrFB',3])
})

#MF - same category
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFS','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propMS']*prop[yr,'propFS'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFs','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrFS',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxFs',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrFS',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propMI']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixFi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMierrFI',3] 
     - simVar[simVar$Year==x & simVar$category=='errMIxFi',3] 
     + simVar[simVar$Year==x & simVar$category=='errMIerrFI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propMB']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMbxFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMberrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMBxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMBerrFB',3])
})

#MF - different categories
#M then F
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propMS']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrFI',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxFi',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrFI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propMS']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrFB',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propMI']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMierrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMIxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMIerrFB',3])
})

#F then M
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propFS']*prop[yr,'propMI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxMi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFserrMI',3] 
     - simVar[simVar$Year==x & simVar$category=='errFSxMi',3] 
     + simVar[simVar$Year==x & simVar$category=='errFSerrMI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propFS']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxMb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFserrMB',3] 
     - simVar[simVar$Year==x & simVar$category=='errFSxMb',3] 
     + simVar[simVar$Year==x & simVar$category=='errFSerrMB',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(1/2)*(1/2)*(prop[yr,'propFI']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFixMb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFierrMB',3] 
     - simVar[simVar$Year==x & simVar$category=='errFIxMb',3] 
     + simVar[simVar$Year==x & simVar$category=='errFIerrMB',3])
})


# Zero SIs. Can't zero FIFB or FSFB as females do contribute to daughters on autosomes
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'MSMI',]
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'FSFI',]
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'MSFI',]
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'FSMI',]
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'MIFI',]
alleleFreqVarAvg_A <-alleleFreqVarAvg[alleleFreqVarAvg$Category != 'MSMI' 
                                    & alleleFreqVarAvg$Category != 'FSFI'
                                    & alleleFreqVarAvg$Category != 'MSFI'
                                    & alleleFreqVarAvg$Category != 'FSMI'
                                  #  & alleleFreqVarAvg$Category != 'MIFI'
                                    ,]

####Z####

ped<-read.table('FSJpedgeno_Zsexlinked.ped',header=FALSE,sep=' ',stringsAsFactors=FALSE)
pedinfo <- ped[,1:5]
colnames(pedinfo) <- c("Family", "USFWS", "Dad", "Mom", "Sex")
indivlist <- merge(simindivFIXmin2obs[,1:6],pedinfo[,c(2,5)],by='USFWS')
indivlist <- indivlist[order(indivlist$Year),]

#samplevar
load("sampleVarZ_21Jan2021_rose.rdata")

#simvar
load("simVarZ_23Jan2021_rose.rdata")


## Number of all & genotyped indivs; proportion of indivs in each category

counts<-ddply(indivlist,.(Year,category,Sex),summarize,genotyped=sum(genotyped=='Y'),
              total=length(category))
counts[counts$Sex==1,'genotyped']<-2*counts$genotyped[counts$Sex==1]
counts[counts$Sex==1,'total']<-2*counts$total[counts$Sex==1]
countsAll<-ddply(indivlist,.(Year,Sex),summarize,genotyped=sum(genotyped=='Y'),
                 total=length(category))
countsAll[countsAll$Sex==1,'genotyped']<-2*countsAll$genotyped[countsAll$Sex==1]
countsAll[countsAll$Sex==1,'total']<-2*countsAll$total[countsAll$Sex==1]

propMS<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='survivor' & counts$Sex==1,'total']/countsAll[countsAll$Year==x & countsAll$Sex==1,'total'])
propFS<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='survivor' & counts$Sex==2,'total']/countsAll[countsAll$Year==x & countsAll$Sex==2,'total'])
propMI<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='immigrant' & counts$Sex==1,'total']/countsAll[countsAll$Year==x & countsAll$Sex==1,'total'])
propFI<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='immigrant' & counts$Sex==2,'total']/countsAll[countsAll$Year==x & countsAll$Sex==2,'total'])
propMB<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='nestling' & counts$Sex==1,'total']/countsAll[countsAll$Year==x & countsAll$Sex==1,'total'])
propFB<-laply(c(2000:2013), function(x) counts[counts$Year==x & counts$category=='nestling' & counts$Sex==2,'total']/countsAll[countsAll$Year==x & countsAll$Sex==2,'total'])
# I think this is proportion of indivs that are in each category and also genotyped? 
# so it is out of total # indivs for that year
prop<-data.frame(propMS=propMS,propFS=propFS,propMI=propMI,propFI=propFI,propMB=propMB,propFB=propFB)
rownames(prop)<-c(2000:2013)
#rowSums(prop)
#Doesn't add to 1 but I think this is because of ungenotyped indivs


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
  (((2/3)*prop[yr,'propMS'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMs-xMt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errMS-errMT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pMspMterrMSerrMT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/3)*prop[yr,'propFS'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFs-xFt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errFS-errFT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pFspFterrFSerrFT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((2/3)*prop[yr,'propMI'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMi-xMt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errMI-errMT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pMipMterrMIerrMT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/3)*prop[yr,'propFI'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFi-xFt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errFI-errFT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pFipFterrFIerrFT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((2/3)*prop[yr,'propMB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xMb-xMt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errMB-errMT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pMbpMterrMBerrMT',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FBsq','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  (((1/3)*prop[yr,'propFB'])^2)*(sampleVar[sampleVar$Year==x & sampleVar$Category=='xFb-xFt','avg'] 
                                 - simVar[simVar$Year==x & simVar$category=='errFB-errFT',3] 
                                 - 2*simVar[simVar$Year==x & simVar$category=='pFbpFterrFBerrFT',3])
})

#MM
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((2/3)^2)*(prop[yr,'propMS']*prop[yr,'propMI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxMi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrMI',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxMi',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrMI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((2/3)^2)*(prop[yr,'propMS']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxMb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrMB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxMb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrMB',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((2/3)^2)*(prop[yr,'propMI']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixMb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMierrMB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMIxMb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMIerrMB',3])
})

#FF
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/3)^2)*(prop[yr,'propFS']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxFi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFserrFI',3] 
     - simVar[simVar$Year==x & simVar$category=='errFSxFi',3] 
     + simVar[simVar$Year==x & simVar$category=='errFSerrFI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/3)^2)*(prop[yr,'propFS']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFserrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errFSxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errFSerrFB',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*((1/3)^2)*(prop[yr,'propFI']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFixFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFierrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errFIxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errFIerrFB',3])
})

#MF - same category
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFS','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propMS']*prop[yr,'propFS'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFs','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrFS',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxFs',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrFS',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propMI']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixFi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMierrFI',3] 
     - simVar[simVar$Year==x & simVar$category=='errMIxFi',3] 
     + simVar[simVar$Year==x & simVar$category=='errMIerrFI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MBFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propMB']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMbxFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMberrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMBxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMBerrFB',3])
})

#MF - different categories
#M then F
alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propMS']*prop[yr,'propFI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrFI',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxFi',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrFI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MSFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propMS']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMsxFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMserrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMSxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMSerrFB',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='MIFB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propMI']*prop[yr,'propFB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xMixFb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xMierrFB',3] 
     - simVar[simVar$Year==x & simVar$category=='errMIxFb',3] 
     + simVar[simVar$Year==x & simVar$category=='errMIerrFB',3])
})

#F then M
alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMI','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propFS']*prop[yr,'propMI'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxMi','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFserrMI',3] 
     - simVar[simVar$Year==x & simVar$category=='errFSxMi',3] 
     + simVar[simVar$Year==x & simVar$category=='errFSerrMI',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FSMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propFS']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFsxMb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFserrMB',3] 
     - simVar[simVar$Year==x & simVar$category=='errFSxMb',3] 
     + simVar[simVar$Year==x & simVar$category=='errFSerrMB',3])
})

alleleFreqVarAvg[alleleFreqVarAvg$Category=='FIMB','avg']<-laply(c(2000:2013), function(x) {
  yr<-as.character(x)
  2*(2/3)*(1/3)*(prop[yr,'propFI']*prop[yr,'propMB'])*
    (sampleVar[sampleVar$Year==x & sampleVar$Category=='xFixMb','avg'] 
     - simVar[simVar$Year==x & simVar$category=='xFierrMB',3] 
     - simVar[simVar$Year==x & simVar$category=='errFIxMb',3] 
     + simVar[simVar$Year==x & simVar$category=='errFIerrMB',3])
})


# Zero SIs, FIFB, and FSFB
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'MSMI',]
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'FSFI',]
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'MSFI',]
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'FSMI',]
#alleleFreqVarAvg[alleleFreqVarAvg$Category == 'MIFI',]
alleleFreqVarAvg_Z <-alleleFreqVarAvg[alleleFreqVarAvg$Category != 'MSMI' 
                                    & alleleFreqVarAvg$Category != 'FSFI'
                                    & alleleFreqVarAvg$Category != 'MSFI'
                                    & alleleFreqVarAvg$Category != 'FSMI'
                                  #  & alleleFreqVarAvg$Category != 'MIFI'
                                    & alleleFreqVarAvg$Category != 'FIFB'
                                    & alleleFreqVarAvg$Category != 'FSFB',]

####combine A and Z####

alleleFreqVarAvg_Z$sum<-laply(alleleFreqVarAvg_Z$Year,function(x)   sum(alleleFreqVarAvg_Z[alleleFreqVarAvg_Z$Year==x,'avg']))
alleleFreqVarAvg_Z$prop<-alleleFreqVarAvg_Z$avg/alleleFreqVarAvg_Z$sum

alleleFreqVarAvg_A$sum<-laply(alleleFreqVarAvg_A$Year,function(x)   sum(alleleFreqVarAvg_A[alleleFreqVarAvg_A$Year==x,'avg']))
alleleFreqVarAvg_A$prop<-alleleFreqVarAvg_A$avg/alleleFreqVarAvg_A$sum

AZ_AFVA<-rbind.data.frame(
  cbind.data.frame(alleleFreqVarAvg_A,"chrom"="A"),
  cbind.data.frame(alleleFreqVarAvg_Z,"chrom"="Z")
)


####label supercategories####
AZ_AFVA$faccat <- factor(AZ_AFVA$Category,
                         levels = 
                           
                           c("MSsq", "MIsq","MBsq",
                             "MSMB" , "MIMB"  ,"FBsq",
                             "FSsq"  , "FIsq",
                             "FSFB"   , "FIFB"    , 
                             "MSFS"   ,   "MBFB"  ,     "MIFI"   , 
                             "MSFB"    ,   "MIFB"   ,   "FSMB"   ,   "FIMB"     
                           ))

AZ_AFVA$Supercategory<-factor(ifelse(AZ_AFVA$Category=='MIFI' | AZ_AFVA$Category=='MIsq' |  AZ_AFVA$Category=='FIsq','Immigrant',
                                     ifelse(AZ_AFVA$Category=='MIMB' | AZ_AFVA$Category=='FIFB' | AZ_AFVA$Category=='MIFB' |  AZ_AFVA$Category=='FIMB', 'Cov(I,B)',
                                            ifelse( AZ_AFVA$Category=='MBFB' | AZ_AFVA$Category=='MBsq' | AZ_AFVA$Category=='FBsq', 'Birth',
                                                    ifelse(AZ_AFVA$Category=='MSMB' | AZ_AFVA$Category=='FSFB' | AZ_AFVA$Category=='MSFB' | AZ_AFVA$Category=='FSMB', 'Cov(S,B)',
                                                           ifelse(AZ_AFVA$Category=='MSsq' |AZ_AFVA$Category=='FSsq' | AZ_AFVA$Category=='MSFS', 'Survivor', NA))))),
                              levels = c('Immigrant', 'Cov(I,B)', 'Birth','Cov(S,B)','Survivor'))

AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("FSsq","FIsq","FBsq","FIFB","FSFB")] <- "Female"
AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("MSsq","MIsq","MBsq","MIMB","MSMB")] <- "Male"
AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("FIMB","FSMB","MIFB","MSFB","MBFB","MSFS","MIFI")] <- "Cov(F,M)"

AZ_AFVA$Category2 <- AZ_AFVA$Category
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("FSsq","MSsq","MSFS")] <- "S"
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("FIsq","MIsq","MIFI")] <- "I"
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("FBsq","MBsq","MBFB")] <- "B"
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("MIMB","FIFB")] <- "IB"
AZ_AFVA$Category2[AZ_AFVA$Category %in% c("MSMB","FSFB")] <- "SB"

AZ_AFVA$Category3 <- AZ_AFVA$Category2
#AZ_AFVA$Category3[AZ_AFVA$Category %in% c("FIMB","FSMB")] <- "FM"
#AZ_AFVA$Category3[AZ_AFVA$Category %in% c("MIFB","MSFB")] <- "MF"

AZ_AFVA$Category4 <- NA
AZ_AFVA$Category4[AZ_AFVA$Category %in% c("FIMB","FSMB")] <- "FM"
AZ_AFVA$Category4[AZ_AFVA$Category %in% c("MIFB","MSFB")] <- "MF"

library(PNWColors)

varp_title <- expression(paste("Variance ",Delta, "p"))


AZ_AFVA$supercategory <- factor(AZ_AFVA$Supercategory,levels=c("Survivor","Birth","Cov(S,B)","Immigrant","Cov(I,B)"))
AZ_AFVA$sexCat <- factor(AZ_AFVA$SexCat,levels=c("Female" ,"Male"    ,   "Cov(F,M)"))

#fills_to_use <-c(pnw_palettes$Bay[1,3],  pnw_palettes$Winter[1,2], pnw_palettes$Bay[1,2],pnw_palettes$Bay[1,1],pnw_palettes$Winter[1,4],pnw_palettes$Bay[1,4],pnw_palettes$Bay[1,5])
fills_to_use <-c(pnw_palettes$Bay[1,3],pnw_palettes$Winter[1,2],"#DD4124", pnw_palettes$Bay[1,2],pnw_palettes$Bay[1,1],pnw_palettes$Winter[1,4],"#E77A66",pnw_palettes$Bay[1,4],pnw_palettes$Bay[1,5])

#Fig S7
#unique(AZ_AFVA$Category3)
ggplot(data=AZ_AFVA, aes(x=Year, y=prop)) + 
  geom_hline(yintercept = 0,alpha=0.5)+
  
  # geom_ribbon(aes(x = Year, ymin = q5_prop, ymax = q95_prop, alpha = chrom)) +
  #geom_errorbar(aes(ymin=q5_prop, ymax=q95_prop,color=Category3), width=.1,alpha=0.5) +
  
  geom_point(aes(color=Category3),alpha=0.5 ) + 
  geom_line(aes(linetype=chrom,color=Category3)) + 
  guides(color=FALSE,linetype=FALSE) +
  #facet_grid(SexCat ~ Supercategory,scales="free")+
  facet_grid(  supercategory ~sexCat,scales="free")+
  
  theme_bw(base_size = 16) + 
  scale_x_continuous(breaks=c(2000,2005,2010),limits = c(1999,2015))+  
  scale_color_manual(
    #breaks=c("FM","MF"),
    values=fills_to_use)+  
  
  labs(y=varp_title,x="Year",color="Category",linetype="") + 
  theme(strip.background =element_rect(fill="white")) +
  geom_dl(aes(label = Category4,color=Category3), method = list("last.qp",cex = 1,dl.trans(x = x + .2))) +
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank()) 

##
AZ_AFVA$year_adj <- AZ_AFVA$Year 
AZ_AFVA$year_adj[AZ_AFVA$chrom == "A"] <- AZ_AFVA$year_adj[AZ_AFVA$chrom == "A"] + 0.25

ggplot(data=AZ_AFVA, aes(x=year_adj, y=prop)) + 
  geom_hline(yintercept = 0,alpha=0.5)+
  
  geom_line(aes(linetype=chrom,color=Category2),alpha=0.5) + 
  
  geom_point(aes(color=Category2,shape=chrom)) + 
  guides(color=FALSE) +
  facet_grid(  supercategory ~sexCat,scales="free")+
  
  theme_bw(base_size = 16) + 
  scale_x_continuous(breaks=c(2000,2005,2010),limits = c(1999,2015))+  
  scale_color_manual(values=fills_to_use)+  
  scale_shape_manual(values=c(1,16)) +
  # scale_linetype_manual(values=c("solid","dashed")) +
  labs(y=varp_title,x="Year",color="Category",linetype="") + 
  theme(strip.background =element_rect(fill="white")) +
  geom_dl(aes(label = Category4,color=Category2), method = list("last.qp",cex = 1,dl.trans(x = x + .2))) +
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank()) #+ 
coord_cartesian(ylim = c(-0.25, 0.5))
