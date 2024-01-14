
SNPlist <- fread('~/Documents/Github/genomeContent/genome_files/SNPsOfInterest.txt')
SNPlist$SNP <- paste0("SNP",genelist$V1)

annotateMySNPs <- function(rec_cm,marey,window_size=1000000,SNPlist=NA){
  
  #get recombination rate
  rec_windows_cm <- rec_cm[[1]]
  rec_windows_cm$logCMMB <- log(rec_windows_cm$cmmb,2)
  
  #get recombination rate quantiles (log scale is much more normal distribution)
  cmquantiles <- quantile(rec_windows_cm$logCMMB,na.rm=T)
  
  rec_windows_cm$cmmbquantile[rec_windows_cm$logCMMB <= cmquantiles[5]] <- "very high"
  rec_windows_cm$cmmbquantile[rec_windows_cm$logCMMB <= cmquantiles[4]] <- "high"
  rec_windows_cm$cmmbquantile[rec_windows_cm$logCMMB <= cmquantiles[3]] <- "low"
  rec_windows_cm$cmmbquantile[rec_windows_cm$logCMMB <= cmquantiles[2]] <- "very low"
  
  print(
    ggplot(rec_windows_cm,aes(x=cmmb,fill=cmmbquantile)) + geom_histogram() + theme_bw()
  )
  
  #round to join to windows
  SNPmap <- marey[,c("SNPname","homolog","mb")]
  
  SNPmap$mb_r <- round(SNPmap$mb,log(window_size/1000000,10))
  SNP_inWins <- left_join(SNPmap,rec_windows_cm,by=c("homolog"="LG","mb_r"="mb_e"))
  
  #take in SNP annotation, usually from SNPeff
  SNP_ann <- fread('~/Documents/Github/genomeContent/genome_files/FSJ2.SNP.annotation.txt')
  SNP_ann$mb <- SNP_ann$V6 / 1000000
  SNP_ann_win <- left_join(SNP_ann,size_homology,by=c("V5"="scaff")) %>% left_join(SNP_inWins,by=c("homolog"="homolog","mb"="mb"))

  SNP_info <- SNP_ann_win %>% dplyr::select(SNPname,V1,homolog,mb,V2,V4,cmmb,cmmbquantile)
  names(SNP_info) <- c("marker","SNP","chr","mb","SNP_impact","geneName","cmmb","logcmmb_quantile")
  
  if(!is.na(SNPlist)){
    SNP_info <- SNP_info[SNP_info$SNP %in% SNPlist$SNP,]
  }
  
  return(SNP_info)
  
  write.table(SNP_info, file = "~/Documents/Github/genomeContent/genome_files/FSJ2.SNP.info.txt", 
              append = F, quote = FALSE, sep = "\t", 
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
              col.names = TRUE, qmethod = c("escape", "double"))
}

makefauxVCF <- function(marey,vcf='~/FSJfullPedFiltDogFINAL12July2016final.vcf'){
  ##function to convert chipseq vcf (usually from plink first) to something that can be annotated
  
  #import vcf
  realvcf <- fread(marey,skip = "#")
  
  #remove individuals
  realvcf10 <- realvcf[,c(1:10)]
  realvcf10_names <- names(realvcf10) #save names for later
  
  #join mapped markers
  SNPmap <- marey[,c("SNPname","NewScaff","SNPpos")]
  fauxvcf <- left_join(realvcf10,SNPmap,by=c("ID"="SNPname"))
  
  #remove old columns
  fauxvcf_reorder <- fauxvcf[!is.na(NewScaff),c("NewScaff","SNPpos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","1_963-80576")]
  
  #sort by chromosome and position
  fauxvcf_sorted <- fauxvcf_reorder[order(fauxvcf_reorder$NewScaff,fauxvcf_reorder$SNPpos),]
  
  #give back old names
  names(fauxvcf_sorted) <- c(realvcf10_names[-10],"fauxind")
  
  #give "the" individual a het site everywhere
  fauxvcf_sorted$fauxind <- "0/1"
  
  write.table(fauxvcf_sorted, file = "~/FSJfullPedFiltDogFINAL12July2016final.faux.vcf", 
              append = F, quote = FALSE, sep = "\t", 
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
              col.names = TRUE, qmethod = c("escape", "double"))
  
  }
