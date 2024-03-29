#script to generate cohort file for genedropping simulations
#Nancy Chen & Rose Driscoll
#Last updated: 03 October 2021

library(plyr)

#read in input file
#nestling data
indiv<-read.table('working_files/IndivData.txt',header=TRUE,sep=' ',stringsAsFactors=FALSE)
indiv$Indiv <- as.character(indiv$Indiv)
cohort<-indiv[indiv$CoreNestling=='Y' & !is.na(indiv$NatalYear),1:2]
write.table(cohort,file="working_files/coreDemoNestlings.txt",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
