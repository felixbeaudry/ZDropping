## Plotting results of A & Z allele frequency model + simulations to estimate error
##  Felix Beaudry and Rose Driscoll 
## June 12 2021
## note: add correct dates to input file names before running

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
library(foreach)

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

####functions####
bootstrapper <- function(allVar){
  col_sample <- sample(length(allVar) -2, 1000, replace=T) + 2
  
  allVar_s <-  allVar[,col_sample]
  allVar_q <- allVar[,c(1:2)]
  
  allVar_q$q5 <- apply(allVar_s, 1, function(x) quantile(x,.05,na.rm = T))
  allVar_q$q95 <- apply(allVar_s, 1, function(x) quantile(x,.95,na.rm = T))
  allVar_q$se <- apply(X=allVar_s,1,function(x) sd(x)/sqrt(length(x)))
  return(allVar_q)
}

popCounter <- function(sex,indivlist,chrom){
  if(sex==T){
    counts<-ddply(indivlist,.(Year,Category,Sex),summarize,Genotyped=sum(Genotyped=='Y'),
                  total=length(Category))
  }
  if(sex==F){
    counts<-ddply(indivlist,.(Year,Sex),summarize,Genotyped=sum(Genotyped=='Y'),
                  total=length(Category))
  }
  if(chrom=="A"){
    counts['Genotyped']<-2*counts$Genotyped
    counts['total']<-2*counts$total
  }
  if(chrom=="Z"){
    counts[counts$Sex==1,'Genotyped']<-2*counts$Genotyped[counts$Sex==1]
    counts[counts$Sex==1,'total']<-2*counts$total[counts$Sex==1]
  }
  if(chrom=="X"){
    counts[counts$Sex==2,'Genotyped']<-2*counts$Genotyped[counts$Sex==2]
    counts[counts$Sex==2,'total']<-2*counts$total[counts$Sex==2]
  }
  return(counts)
}

propsTabler <- function(indivlistgeno,Year=c(2000:2013),chrom){
  indivlist <- indivlistgeno[,c(1:7)]
  names(indivlist)[c(1:7)] <- c( "Year","Indiv", "Category", "Genotyped", "Mom", "Dad", "Sex")
  
  counts <- popCounter(sex=T,indivlist,chrom)
  countsAll <- popCounter(sex=F,indivlist,chrom)
  
  #proportion of indivs that are in each category and also Genotyped 
  propMS<-laply(Year, function(x) counts[counts$Year==x & counts$Category=='survivor' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
  propFS<-laply(Year, function(x) counts[counts$Year==x & counts$Category=='survivor' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
  propMI<-laply(Year, function(x) counts[counts$Year==x & counts$Category=='immigrant' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
  propFI<-laply(Year, function(x) counts[counts$Year==x & counts$Category=='immigrant' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
  propMB<-laply(Year, function(x) counts[counts$Year==x & counts$Category=='nestling' & counts$Sex==1,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
  propFB<-laply(Year, function(x) counts[counts$Year==x & counts$Category=='nestling' & counts$Sex==2,'total']/(countsAll[countsAll$Year==x & countsAll$Sex==1,'total']+countsAll[countsAll$Year==x & countsAll$Sex==2,'total']))
  
  prop <-data.frame(propMS=propMS,propFS=propFS,propMI=propMI,propFI=propFI,propMB=propMB,propFB=propFB)
  rownames(prop)<-Year
  return(prop)
}

alleleFrequer_sq <- function(Category,propterm,sampleCat,errCat,errCov,core_data,boot_data,prop_data,Year = c(2000:2013)){
  
  avg <- 
    laply(Year, function(x) {
      yr<-as.character(x)
      ((prop_data[yr,propterm])^2)*(core_data[core_data$Year==x & core_data$Category==sampleCat & core_data$type == 'sample','avg'] 
                               - core_data[core_data$Year==x & core_data$Category==errCat & core_data$type == 'sim',3] 
                               - 2*core_data[core_data$Year==x & core_data$Category==errCov & core_data$type == 'sim',3])
    })
  
  quantiles <- 
    foreach(quantile= c('q5','q95'),.combine=cbind) %do% {
      laply(Year, function(x) {
        yr<-as.character(x)
        ((prop_data[yr,propterm])^2)*(boot_data[boot_data$Year==x & boot_data$Category==sampleCat,quantile] 
                                 - boot_data[boot_data$Year==x & boot_data$Category==errCat,quantile] 
                                 - 2*boot_data[boot_data$Year==x & boot_data$Category==errCov,quantile])
      })
      
    }
  
  return(cbind(Category,Year,avg,quantiles))  
  
}


alleleFrequer_cov <- function(Category,prop1,prop2,sampleCat,errCat1,errCat2,errCov,prop_data,core_data,boot_data,Year = c(2000:2013)){
  
  avg <- 
    laply(Year, function(x) {
      yr<-as.character(x)
      
      2*(prop_data[yr,prop1]*prop_data[yr,prop2])*
        (core_data[core_data$Year==x & core_data$Category==sampleCat & core_data$type == 'sample','avg'] 
         - core_data[core_data$Year==x & core_data$Category==errCat1 & core_data$type == 'sim',3] 
         - core_data[core_data$Year==x & core_data$Category==errCat2 & core_data$type == 'sim',3]
         + core_data[core_data$Year==x & core_data$Category==errCov & core_data$type == 'sim',3]
        )
      
    })
  
  quantiles <- 
    foreach(quantile= c('q5','q95'),.combine=cbind) %do% {
      serrCat1 <- paste0('s',errCat1)
      laply(Year, function(x) {
        yr<-as.character(x)
        2*(prop_data[yr,prop1]*prop_data[yr,prop2])*
          (boot_data[boot_data$Year==x & boot_data$Category==sampleCat ,quantile] 
           - boot_data[boot_data$Year==x & boot_data$Category==serrCat1 ,quantile] 
           - boot_data[boot_data$Year==x & boot_data$Category==errCat2,quantile]
           + boot_data[boot_data$Year==x & boot_data$Category==errCov ,quantile]
          )
        
      })
      
    }
  
  return(cbind(Category,Year,avg,quantiles))  
  
}
alleleFrequer_full <- function(categories_table_sq,categories_table_cov,core_data,boot_data,prop_data){
  categories_table <- read.table(categories_table_sq,header = T)
  categories_table_cov <- read.table(categories_table_cov,header = T)
  
  alleleFreqVarAvg_sq <- foreach(n= c(1:nrow(categories_table)),.combine=rbind) %do% {  
    print(categories_table[n,1])
    alleleFrequer_sq(
      Category=categories_table[n,1],
      propterm=categories_table[n,2],
      sampleCat=categories_table[n,3],
      errCat=categories_table[n,4],
      errCov=categories_table[n,5],
      
      core_data=core_data,
      boot_data=boot_data,
      prop_data=prop_data)
  }
  
  alleleFreqVarAvg_cov <- foreach(n= c(1:nrow(categories_table_cov)),.combine=rbind) %do% {  
    print(categories_table_cov[n,1])
    alleleFrequer_cov(Category=categories_table_cov[n,1],
                      prop1=categories_table_cov[n,2],
                      prop2=categories_table_cov[n,3],
                      sampleCat=categories_table_cov[n,4],
                      errCat1=categories_table_cov[n,5],
                      errCat2=categories_table_cov[n,6],
                      errCov=categories_table_cov[n,7],
                      core_data=core_data,
                      boot_data=boot_data,
                      prop_data=prop_data)
    
    
  }
  
  alleleFreqVarAvg <- rbind.data.frame(
    alleleFreqVarAvg_sq,
    alleleFreqVarAvg_cov
  )
  
  names(alleleFreqVarAvg)[c(4,5)] <- c("q5","q95")
  alleleFreqVarAvg[c(2:5)] <- sapply(alleleFreqVarAvg[c(2:5)],as.numeric)
  return(alleleFreqVarAvg)
}  

freq_normalizer <- function(alleleFreqVarAvg){
  # Find sum for each year and calculate proportion for each Category
  alleleFreqVarAvg$sum<-laply(alleleFreqVarAvg$Year,function(x)   sum(alleleFreqVarAvg[alleleFreqVarAvg$Year==x,'avg']))
  
  alleleFreqVarAvg$prop<-alleleFreqVarAvg$avg/alleleFreqVarAvg$sum
  
  alleleFreqVarAvg$q5_prop <- alleleFreqVarAvg$q5 / alleleFreqVarAvg$sum
  alleleFreqVarAvg$q95_prop <- alleleFreqVarAvg$q95 / alleleFreqVarAvg$sum
  return(alleleFreqVarAvg)
}

alleleFrequer_birth <- function(Category,propCat,sampleCat,errCat,core_data,prop_data,Year = c(2000:2013)){
  
  avg <- 
    laply(Year, function(x) {
      yr<-as.character(x)
      
      2*(prop_data[yr,propCat])*
        (core_data[core_data$Year==x & core_data$Category==sampleCat & core_data$type == 'sample',3] 
         - core_data[core_data$Year==x & core_data$Category==errCat & core_data$type == 'sim',3] 
        )
      
    })
  return(cbind(Category,Year,avg))  
} 

mendAndFamer <- function(birthCatFile,core_data,prop_data){
  categories_table_birth <- read.table(birthCatFile,header = T)
  
  newBsq <- foreach(n= c(1:nrow(categories_table_birth)),.combine=rbind) %do% {  
    print(categories_table_birth[n,1])
    alleleFrequer_birth(Category=categories_table_birth[n,1],
                        propCat=categories_table_birth[n,2],
                        sampleCat=categories_table_birth[n,3],
                        errCat=categories_table_birth[n,4],
                        core_data=core_data,
                        prop_data=prop_data)
    
  }
  newBsq <- as.data.frame(newBsq)
  newBsq[c(2:3)] <- sapply(newBsq[c(2:3)],as.numeric)
  return(newBsq)
}

OBSvEXPer <- function(core_data,alleleFreqVarAvg,chrom){
  obsDiff <- laply(c(2000:2013), function(x) {
    yr<-as.character(x)
    cbind(yr,(core_data[core_data$Year==x & core_data$Category=='xt1-xt' & core_data$type == 'sample','avg'] 
              - core_data[core_data$Year==x & core_data$Category=='errt1-errt' & core_data$type == 'sim',3]
              - 2*core_data[core_data$Year==x & core_data$Category=='pt1pterrt1errT' & core_data$type == 'sim',3]
    ))
  })
  
  obsDiff <- as.data.frame(obsDiff)
  obsDiff$yr <- as.numeric(obsDiff$yr)
  obsDiff$V2 <- as.numeric(obsDiff$V2)
  
  expDiff <- unique(alleleFreqVarAvg[,c(2,6)])
  
  EO <- left_join(expDiff,obsDiff,by=c("Year"="yr"))
  
  print(
    ggplot(EO,aes(x=sum,y=V2)) + 
      geom_abline(slope=1,color="grey",linetype="dashed") +
      geom_point() + theme_bw() + 
      labs(title=chrom,x=expression(paste("Predicted Change ",Sigma)),y=expression(paste("Expected Change (", p["t"],"-", p["t-1"], ")"))) + 
      theme(aspect.ratio = 1)
  )
  print(
    cor(EO$V2, EO$sum)
  )
  EO_lm <- lm(V2~ sum,data=EO)
  print(
    summary(EO_lm)
  )
}


####start calculations for autosomal loci####

load(file='../working_files/intermediate_files/indivlistgeno_A.rdata')

## Number of all & Genotyped indivs; proportion of indivs in each category
prop_A <- propsTabler(indivlistgeno_A,chrom = "A")

#cat sample and sim tables
load("sampleVarA_pruned.rdata") #samplevar
load("simVarA_pruned.rdata") #simvar

names(simVar) <- names(sampleVar)
simVar$type <- "sim"
sampleVar$type <- "sample"

core_data_A <- rbind(sampleVar,simVar) 

#bootstrap from windows
load("allVar_boot_A_w3.4mb_pruned.rdata") #bootstrap
allVar_q_A <- bootstrapper(allVar)

#calculate final terms for each year, according to formula: (roughly) fraction of population * (change in frequency - error)
#also calculate bootstrap quantiles 
alleleFreqVarAvg_A_raw <- alleleFrequer_full(categories_table_sq='categories_table.txt',
                                             categories_table_cov='categories_table_cov.txt',
                                             core_data=core_data_A,
                                             boot_data=allVar_q_A,
                                             prop_data=prop_A)


  
#remove irrelevant covariance terms
alleleFreqVarAvg_A <-alleleFreqVarAvg_A_raw[
                                      alleleFreqVarAvg_A_raw$Category != 'MSMI' 
                                       & alleleFreqVarAvg_A_raw$Category != 'FSFI'
                                       & alleleFreqVarAvg_A_raw$Category != 'MSFI'
                                       & alleleFreqVarAvg_A_raw$Category != 'FSMI'
                                       ,]
  
#normalize terms#
alleleFreqVarAvg_A <-  freq_normalizer(alleleFreqVarAvg_A)

#Mend and Fam
newBsq_A <- mendAndFamer(birthCatFile='categories_table_birth.txt',core_data=core_data_A,prop_data=prop_A)

OBSvEXPer(core_data = core_data_A,alleleFreqVarAvg=alleleFreqVarAvg_A,chrom="A")

####start Z####

#indivlist
load(file='../working_files/intermediate_files/indivlistgeno_Z.rdata')

prop_Z <- propsTabler(indivlistgeno_Z,chrom="Z")


load("sampleVarZ_pruned.rdata") #samplevar
load("simVarZ_pruned.rdata") #simvar

names(simVar) <- names(sampleVar)
simVar$type <- "sim"
sampleVar$type <- "sample"

core_data_Z <- rbind(sampleVar,simVar) 


load("allVar_boot_Z_w3.4mb_pruned.rdata") #bootstrap

allVar_q_Z <- bootstrapper(allVar)

alleleFreqVarAvg_Z_raw <- alleleFrequer_full(categories_table_sq='categories_table.txt',
                                             categories_table_cov='categories_table_cov.txt',
                                             core_data=core_data_Z,
                                             boot_data=allVar_q_Z,
                                             prop_data=prop_Z)



alleleFreqVarAvg_Z<-alleleFreqVarAvg_Z_raw[
                              
                                      alleleFreqVarAvg_Z_raw$Category != 'MSMI' 
                                    & alleleFreqVarAvg_Z_raw$Category != 'FSFI'
                                    & alleleFreqVarAvg_Z_raw$Category != 'MSFI'
                                    & alleleFreqVarAvg_Z_raw$Category != 'FSMI'
                                    
                                    # Remove immigrant covariance and mom-daughter Z inheritance
                                    
                                    & alleleFreqVarAvg_Z_raw$Category != 'FIFB'
                                    & alleleFreqVarAvg_Z_raw$Category != 'FSFB',]


alleleFreqVarAvg_Z <-  freq_normalizer(alleleFreqVarAvg_Z)

newBsq_Z <- mendAndFamer(birthCatFile='categories_table_birth.txt',
                         core_data=core_data_Z,prop_data=prop_Z)


OBSvEXPer(core_data = core_data_Z,alleleFreqVarAvg=alleleFreqVarAvg_Z,chrom="Z")


####combine A and Z####

#stack
AZ_AFVA<-rbind.data.frame(
  cbind.data.frame(alleleFreqVarAvg_A,"chrom"="A"),
  cbind.data.frame(alleleFreqVarAvg_Z,"chrom"="Z")
)

write.table(AZ_AFVA, file = "AZ_AFVA_prune.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

####label supercategories for stack####
#detach(package:plyr)

#plotting parameters
Avarp_title <- expression(paste("Autosomal Variance ",Delta, "p"))
Zvarp_title <- expression(paste("Z Variance ",Delta, "p"))
varp_title <- expression(atop("Proportional Contribution",paste("to ",Delta, "p variance")))
varp_title2 <- expression(paste("Proportional Contribution to ",Delta, "p variance"))

fills_to_use <-c(pnw_palettes$Bay[1,3],pnw_palettes$Winter[1,2],"#DD4124", pnw_palettes$Bay[1,2],pnw_palettes$Bay[1,1],pnw_palettes$Winter[1,4],"#E77A66",pnw_palettes$Bay[1,4],pnw_palettes$Bay[1,5])

####Plot Fig. 5####

####categories and sexes####

#set sex label & factor
AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("FSsq","FIsq","FBsq","FIFB","FSFB")] <- "Female"
AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("MSsq","MIsq","MBsq","MIMB","MSMB")] <- "Male"
AZ_AFVA$SexCat[AZ_AFVA$Category %in% c("FIMB","FSMB","MIFB","MSFB","MBFB","MSFS","MIFI")] <- "Cov(F,M)"
AZ_AFVA$sexCat <- factor(AZ_AFVA$SexCat,levels=c("Female" ,"Male"    ,   "Cov(F,M)"))


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
AZ_AFVA$year_adj <- AZ_AFVA$Year 
AZ_AFVA$year_adj[AZ_AFVA$chrom == "A"] <- AZ_AFVA$year_adj[AZ_AFVA$chrom == "A"] + 0.25


pdf(paste("fig_5.prune.pdf",sep=''),width=5.5,height=4)
base_font_size = 8

plot_grid(

ggplot(data=AZ_AFVA %>% filter(supercategory == "Survivor"), aes(x=year_adj, y=prop)) + 
  geom_hline(yintercept = 0,alpha=0.5,size=0.25)+
  
  geom_errorbar(aes(ymin=q5_prop, ymax=q95_prop,color=Category2,linetype=chrom), width=.1,alpha=0.5,size=0.15) +
  #  geom_line(aes(linetype=chrom,color=Category2),alpha=0.5) + 
  
  geom_point(aes(color=Category2,shape=chrom),size=0.75) + 
  facet_grid(  supercategory ~sexCat,scales="free")+
  
  theme_bw(base_size = base_font_size) + 
  scale_x_continuous(breaks=c(2000,2005,2010),limits = c(1999,2015))+  
  scale_color_manual(values=fills_to_use[1])+  
  scale_shape_manual(values=c(16,1)) +
  # scale_linetype_manual(values=c("solid","dashed")) +
#  labs(y=varp_title2,x="Year",color="Category",linetype="") + 
  theme(strip.background =element_rect(fill="white")) +
  geom_dl(aes(label = Category4,color=Category2), method = list("last.qp",cex = 0.5,dl.trans(x = x + .2))) +
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  coord_cartesian(ylim = c(-0.3, 0.5))  +
  
  guides(color=FALSE,shape=FALSE,linetype=FALSE) +  labs(y="",x="Year",color="Category",linetype="") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )  + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank())
,
ggplot(data=AZ_AFVA %>% filter(supercategory == "Cov(S,B)"), aes(x=year_adj, y=prop)) + 
  geom_hline(yintercept = 0,alpha=0.5,size=0.25)+
  
  geom_errorbar(aes(ymin=q5_prop, ymax=q95_prop,color=Category2,linetype=chrom), width=.1,alpha=0.5,size=0.15) +
  #  geom_line(aes(linetype=chrom,color=Category2),alpha=0.5) + 
  
  geom_point(aes(color=Category2,shape=chrom),size=0.75) + 
  facet_grid(  supercategory ~sexCat,scales="free")+
  
  theme_bw(base_size = base_font_size) + 
  scale_x_continuous(breaks=c(2000,2005,2010),limits = c(1999,2015))+  
  scale_color_manual(values=fills_to_use[c(2:4)])+  
  scale_shape_manual(values=c(16,1)) +
  # scale_linetype_manual(values=c("solid","dashed")) +
  theme(strip.background =element_rect(fill="white")) +
  geom_dl(aes(label = Category4,color=Category2), method = list("last.qp",cex = 0.5,dl.trans(x = x + .2))) +
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  coord_cartesian(ylim = c(-0.2, 0.4)) +
  
  guides(color=FALSE,shape=FALSE,linetype=FALSE) +  labs(y="",x="Year",color="Category",linetype="") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )  + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank())
,
ggplot(data=AZ_AFVA %>% filter(supercategory ==  "Birth"), aes(x=year_adj, y=prop)) + 
  geom_hline(yintercept = 0,alpha=0.5,size=0.25)+
  
  geom_errorbar(aes(ymin=q5_prop, ymax=q95_prop,color=Category2,linetype=chrom), width=.1,alpha=0.5,size=0.15) +
  #  geom_line(aes(linetype=chrom,color=Category2),alpha=0.5) + 
  
  geom_point(aes(color=Category2,shape=chrom),size=0.75) + 
  facet_grid(  supercategory ~sexCat,scales="free")+
  
  theme_bw(base_size = base_font_size) + 
  scale_x_continuous(breaks=c(2000,2005,2010),limits = c(1999,2015))+  
  scale_color_manual(values=fills_to_use[5])+  
  scale_shape_manual(values=c(16,1)) +
  theme(strip.background =element_rect(fill="white")) +
  geom_dl(aes(label = Category4,color=Category2), method = list("last.qp",cex = 0.5,dl.trans(x = x + .2))) +
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  
  coord_cartesian(ylim = c(-0.1, 0.5)) +
  
  guides(color=FALSE,shape=FALSE,linetype=FALSE) +  labs(y="",x="Year",color="Category",linetype="") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )  + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank())
,

ggplot(data=AZ_AFVA %>% filter(supercategory ==  "Cov(I,B)"), aes(x=year_adj, y=prop)) + 
  geom_hline(yintercept = 0,alpha=0.5,size=0.25)+
  
  geom_errorbar(aes(ymin=q5_prop, ymax=q95_prop,color=Category2,linetype=chrom), width=.1,alpha=0.5,size=0.15) +
  #  geom_line(aes(linetype=chrom,color=Category2),alpha=0.5) + 
  
  geom_point(aes(color=Category2,shape=chrom),size=0.75) + 

  facet_grid(  supercategory ~sexCat,scales="free")+
  
  theme_bw(base_size = base_font_size) + 
  scale_x_continuous(breaks=c(2000,2005,2010),limits = c(1999,2015))+  
  scale_color_manual(values=fills_to_use[c(6:8)])+  
  scale_shape_manual(values=c(16,1)) +
  # scale_linetype_manual(values=c("solid","dashed")) +

  theme(strip.background =element_rect(fill="white")) +
  geom_dl(aes(label = Category4,color=Category2), method = list("last.qp",cex = 0.5,dl.trans(x = x + .2))) +
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  
 coord_cartesian(ylim = c(-0.1, 0.2))  +
  
  guides(color=FALSE,shape=FALSE,linetype=FALSE) +  labs(y="",x="Year",color="Category",linetype="") + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )  + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank())
,
ggplot(data=AZ_AFVA %>% filter(supercategory ==  "Immigrant"), aes(x=year_adj, y=prop)) + 
  geom_hline(yintercept = 0,alpha=0.5,size=0.25)+
  
  geom_errorbar(aes(ymin=q5_prop, ymax=q95_prop,color=Category2,linetype=chrom), width=.1,alpha=0.5,size=0.15) +
  #  geom_line(aes(linetype=chrom,color=Category2),alpha=0.5) + 
  
  geom_point(aes(color=Category2,shape=chrom),size=0.75) + 
  
  facet_grid(  supercategory ~sexCat,scales="free")+
  
  theme_bw(base_size = base_font_size) + 
  scale_x_continuous(breaks=c(2000,2005,2010),limits = c(1999,2015))+  
  scale_color_manual(values=fills_to_use[9])+  
  scale_shape_manual(values=c(16,1)) +
  # scale_linetype_manual(values=c("solid","dashed")) +
  theme(strip.background =element_rect(fill="white")) +
  geom_dl(aes(label = Category4,color=Category2), method = list("last.qp",cex = 0.5,dl.trans(x = x + .2))) +
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  
  coord_cartesian(ylim = c(-0.01, 0.08)) +

guides(color=FALSE,shape=FALSE,linetype=FALSE) +  labs(y="",x="Year",color="Category",linetype="") + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
)  + theme(
         axis.text.x = element_blank(),
         axis.title.x = element_blank())


,  labels = NA, ncol = 1, align = 'v',axis='rl')





dev.off()


min(AZ_AFVA$prop[AZ_AFVA$supercategory == "Immigrant"])
max(AZ_AFVA$prop[AZ_AFVA$supercategory == "Immigrant"])

min(AZ_AFVA$prop[AZ_AFVA$supercategory == "Birth"])
max(AZ_AFVA$prop[AZ_AFVA$supercategory == "Birth"])

min(AZ_AFVA$prop[AZ_AFVA$supercategory == "Survivor"])
max(AZ_AFVA$prop[AZ_AFVA$supercategory == "Survivor"])


###cat and sex specific plots

AZ_AFVA_sex <- 
  AZ_AFVA %>% group_by(chrom,Year) %>% 
  mutate(sex_sum=sum(avg)) %>% 
  group_by(chrom,Year,SexCat) %>% 
  mutate(sexCat_sum=sum(avg)/sex_sum) 

AZ_AFVA_sex_only <- unique(AZ_AFVA_sex %>% dplyr::select(Year ,chrom ,sexCat , sex_sum ,sexCat_sum)) #[,c(1,10,12:14)])
AZ_AFVA_sex_only$SexCat = factor(AZ_AFVA_sex_only$sexCat, levels=c('Cov(F,M)','Female','Male'))

AZ_AFVA_cat <- 
  AZ_AFVA %>% group_by(chrom,Year) %>% 
  mutate(cat_sum=sum(avg)) %>% 
  group_by(chrom,Year,supercategory) %>% 
  mutate(catCat_sum=sum(avg)/cat_sum) 


AZ_AFVA_cat_only <- unique(AZ_AFVA_cat %>% dplyr::select( Year, chrom ,supercategory,  cat_sum, catCat_sum)) #[,c(1,10,15,18:19)])
AZ_AFVA_cat_only$supercategory = factor(AZ_AFVA_cat_only$supercategory, levels=c('Cov(I,B)','Immigrant','Cov(S,B)','Birth','Survivor'))

fills_to_use_cat <-c(pnw_palettes$Bay[1,1],pnw_palettes$Bay[1,2],pnw_palettes$Bay[1,5], pnw_palettes$Bay[1,3],pnw_palettes$Bay[1,4])

min(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Male" & AZ_AFVA_sex_only$chrom == "A"])
max(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Male" & AZ_AFVA_sex_only$chrom == "A"])

min(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Female" & AZ_AFVA_sex_only$chrom == "A"])
max(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Female" & AZ_AFVA_sex_only$chrom == "A"])

min(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Male" & AZ_AFVA_sex_only$chrom == "Z"])
max(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Male" & AZ_AFVA_sex_only$chrom == "Z"])

min(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Female" & AZ_AFVA_sex_only$chrom == "Z"])
max(AZ_AFVA_sex_only$sexCat_sum[AZ_AFVA_sex_only$sexCat == "Female" & AZ_AFVA_sex_only$chrom == "Z"])

####plot fig. 4B####
pdf(paste("fig_4BC.prune.pdf",sep=''),width=4.5,height=3)

plot_grid(
#sex_stack_plot <- 
  ggplot() + 
  geom_hline(yintercept = 0,alpha=0.5)+
  geom_hline(yintercept = 1,alpha=0.5)+
  
  geom_bar(data=AZ_AFVA_sex_only, aes(x=Year, y=sexCat_sum,fill=SexCat,group=SexCat), stat='identity',alpha=0.75) +
  
  facet_grid(chrom~.,scales="free",space="free") + 
  labs(y=varp_title,x="Year") +
  theme_bw(base_size = base_font_size) + 
  theme(strip.background =element_rect(fill="white")) +
  scale_fill_manual("Category",values=c("mediumpurple","indianred1","cornflowerblue")) + 
  guides(color=FALSE)  + 
  theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_y_continuous(breaks=c(0,0.5,1))
,
#cat_stack_plot <- 
  ggplot() + 
  geom_hline(yintercept = 0,alpha=0.5)+
  
  geom_hline(yintercept = 1,alpha=0.5)+
  
  geom_bar(data=AZ_AFVA_cat_only, aes(x=Year, y=catCat_sum,fill=supercategory,group=supercategory), stat='identity',alpha=0.75) +
  
  facet_grid(chrom~.,scales="free",space="free") + 
  labs(y=varp_title,x="Year") +
  theme_bw(base_size = base_font_size) + 
  theme(strip.background =element_rect(fill="white")) +
  scale_fill_manual("Category",values=fills_to_use_cat) + 
  guides(color=FALSE)  + 
  theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_y_continuous(breaks=c(0,0.5,1))


,  labels = c(' ', ' '), ncol = 1, align = 'v',axis='rl')

dev.off()

##
rep_sum_A <- AZ_AFVA_cat_only[AZ_AFVA_cat_only$supercategory %in% c("Survivor","Birth","Cov(S,B)") & AZ_AFVA_cat_only$chrom == "A",] %>% 
  group_by(chrom,Year) %>% 
  summarize(rep_sum=sum( catCat_sum)) 

rep_sum_Z <- AZ_AFVA_cat_only[AZ_AFVA_cat_only$supercategory %in% c("Survivor","Birth","Cov(S,B)") & AZ_AFVA_cat_only$chrom == "Z",] %>% 
  group_by(chrom,Year) %>% 
  summarize(rep_sum=sum( catCat_sum)) 

####label supercategories for merge####
alleleFreqVarAvg1_AZ <- left_join(alleleFreqVarAvg_A,alleleFreqVarAvg_Z,by=c("Year"="Year","Category"="Category"))

names(alleleFreqVarAvg1_AZ) <- 
c("Category" ,"Year"   ,    "avg.A"   ,   "q5.A"   ,    "q95.A"    ,  "sum.A"  ,   
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


 pdf(paste("fig_S7.pdf",sep=''),width=6,height=4)
 
 ggplot(alleleFreqVarAvg1_AZ[alleleFreqVarAvg1_AZ$FullCat %ni% c("Cov(Fs,Fb)","Cov(Fi,Fb)"),],aes(x=prop.A,y=prop.Z)) +   
  geom_abline(intercept = 0,slope=(4/3),alpha=0.5,linetype="dotted",size=0.25,color="darkgreen") +
  geom_abline(intercept = 0,slope=1,alpha=0.5,size=0.25,color="darkgreen") +
  geom_abline(intercept = 0,slope=(2/3),alpha=0.5,linetype="dashed",size=0.25,color="darkgreen") +
  
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

####OBS EXP####




