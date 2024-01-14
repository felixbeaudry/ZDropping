#R script to build annotation in windows, with weighting based on BUSCO set
#F. Beaudry as of feb 6, 22

library(tidyverse)
library(data.table)
require(R.utils) #for opening gzipped files

library(foreach)
library(doParallel)

source("genomeContent_functions.R")

options(scipen=999)

cores=detectCores() #uncomment these two lines if you want to use more than 4 cores
cl <- makeCluster(round(cores[1]/2)) #not to overload your computer
#cl <- makeCluster(4) #not to overload your computer
registerDoParallel(cl)

#setwd('~/Documents/Github/genomeContent/')

coding_window_size = 1000
windowSizelog = log(0.0001,10) 
braker_score_cutoff = 0.85
aves_busco_N = 8338

today <- "220226"

#windowSizelog number is for the size of the window to merge identical busco genes, in mb so 0.0001 = 100bp 
#window sized determined by answering the question 'what size leads to accurate merging of tandem duplicates?' manually

####functions####


####file intake####
homology_convert <- fread('genome_files/FSJ2.convert.txt')
sizes <- fread('genome_files/FSJ2.scaffSizes.txt')
homology_sizes <- left_join(sizes,homology_convert,by=c("V1"="V2"))
names(homology_sizes) <- c("scaf","chrom_length","homolog")

genome_busco <- getGenomeBUSCO(genomeBUSCOFile='genome_files/annotation_files/FSJ2.genome_busco.tsv')

stringtie <- reformatGTF(predictionType="stringtie",
                                  inputFile= 'genome_files/annotation_files/FSJ2.stringtieMerged.gtf',
                                  buscoFile="genome_files/annotation_files/FSJ2.stringtie_busco.tsv")

braker1 <- reformatGTF(predictionType="braker1",
                       inputFile='genome_files/annotation_files/FSJ2.BRAKER1.gtf.gz',
                       buscoFile="genome_files/annotation_files/FSJ2.BRAKER1_busco.tsv")


braker2 <- reformatGTF(predictionType="braker2",
                        inputFile='genome_files/annotation_files/FSJ2.BRAKER2.gtf.gz',
                       buscoFile="genome_files/annotation_files/FSJ2.BRAKER2_busco.tsv")

genes <- rbind.data.frame(braker1,braker2,stringtie)
genes$gene_length <- genes$end_pos - genes$start_pos
genes <- genes[order(genes$scaffold,genes$start_pos),]

####windows####

coding_win_loop_chr <- 
  foreach(chrom=homology_sizes$scaf[!is.na(homology_sizes$homolog)][c(22)],.combine=rbind) %dopar% {
    #chrom = "ScYP8k369HRSCAF175ch21"
    library(tidyverse)
    library(foreach)

  win_loop_tmp <- 
    foreach(win=seq(0,coding_window_size*100,coding_window_size),.combine=rbind) %do% {
      
  #  foreach(win=seq(0,homology_sizes$chrom_length[homology_sizes$scaf == chrom],coding_window_size),.combine=rbind) %do% {
      #win=0
        
        #      [1]     [2]     [3]        [4]             [5]             [6]             [7]         [8]             [9]             [10]             [11]             [12]             [13]             [14]
        # c("chrom","win_s","win_e","braker1_gene","braker1_score","braker1_busco","braker1_bp", "braker2_gene","braker2_score","braker2_busco","braker2_bp","stringtie_gene","stringtie_busco","stringtie_bp")
        tmp_frame <- c(chrom,win,(win+coding_window_size),NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
       
        genes_tmp <- genes %>% filter(scaffold == chrom & start_pos < (win+coding_window_size) & end_pos > win )
        
        for(thisOrigin in unique(genes_tmp$origin)){
          
          colPush = 0 
          if(thisOrigin == "braker2" ){colPush = 4}
          if(thisOrigin == "stringtie" ){colPush = 8}
          
          genes_byOrigin <- genes_tmp %>% filter(origin == thisOrigin)
          numberOfGenes <- length(unique(genes_byOrigin$gene_id))
          
          if(numberOfGenes == 1){
            tmp_frame[(4 + colPush )] <- unique(genes_byOrigin$gene_id)
            tmp_frame[(4 + colPush + 1)] <- max(genes_byOrigin$annotation_score)
            tmp_frame[(4 + colPush + 2)] <- unique(genes_byOrigin$busco_id) #busco
            tmp_frame[(4 + colPush + 3)] <- returnLongestTranscript(genes_byOrigin,win,coding_window_size)[1,2] #length
          }
          
          if(numberOfGenes > 1){
            coding_length = 0 
            
            for(thisGene in unique(genes_byOrigin$gene_id)){
              genes_byGene <- genes_byOrigin %>% filter(gene_id == thisGene)
              numberOfTranscripts <- length(unique(genes_byGene$transcript_id))
              
              longestTranscript <- returnLongestTranscript(genes_byGene,win,coding_window_size)
              coding_length = coding_length + longestTranscript$codingLength
              
            }
            
            tmp_frame[(4 + colPush )] <- "overlap"
            #tmp_frame[(4 + colPush + 1)] <- max(genes_byOrigin$annotation_score)
            #tmp_frame[(4 + colPush + 2)] <- unique(genes_byOrigin$busco_id) #busco
            tmp_frame[(4 + colPush + 3)] <- coding_length #length
            
          }
        }
          
        tmp_frame
      }
    
    win_loop_tmp_df <- as.data.frame(win_loop_tmp)
    names(win_loop_tmp_df) <-      c("chrom","win_s","win_e","braker1_gene","braker1_score","braker1_busco","braker1_bp", "braker2_gene","braker2_score","braker2_busco","braker2_bp","stringtie_gene","stringtie_score","stringtie_busco","stringtie_bp")
    
    win_loop_tmp_df
}  
stopCluster(cl)

coding_win_loop_chr[, c(2,3,5,7,9,11,13,15)] <- sapply(coding_win_loop_chr[, c(2,3,5,7,9,11,13,15)], as.numeric)
coding_win_loop_chr[, -c(2,3,5,7,9,11,13,15)] <- sapply(coding_win_loop_chr[, -c(2,3,5,7,9,11,13,15)], as.character)
#coding_win_loop_chr[, c(2,3,5,7,9,11,14)] <- sapply(coding_win_loop_chr[, c(2,3,5,7,9,11,14)], as.numeric)
#coding_win_loop_chr[, -c(2,3,5,7,9,11,14)] <- sapply(coding_win_loop_chr[, -c(2,3,5,7,9,11,14)], as.character)

coding_win_loop_chr[coding_win_loop_chr == "NA"] <- NA

coding_win <- as.data.frame(coding_win_loop_chr)

save(coding_win,file="genome_files/annotation_files/FSJ2.coding_win.rdata")


####gene content####
load(file="genome_files/annotation_files/FSJ2.coding_win.rdata")

coding_win$braker1_gene[coding_win$braker1_score < braker_score_cutoff] <- "bad_exon"
coding_win$braker2_gene[coding_win$braker2_score < braker_score_cutoff] <- "bad_exon"

coding_win$braker1_score[is.na(coding_win$braker1_score)] <- 0
coding_win$braker2_score[is.na(coding_win$braker2_score)] <- 0

coding_win$braker_best <- ifelse(coding_win$braker1_score >= coding_win$braker2_score,coding_win$braker1_gene,coding_win$braker2_gene)
coding_win$braker_best[coding_win$braker_best == "bad_exon"] <- NA

braker_best <- findBestHomolog(coding_win=coding_win,homolog_mode="braker")

coding_win<- left_join(coding_win,braker_best,by=c("braker_best"="worst"))

coding_win$braker_best[!is.na(coding_win$best)] <- coding_win$best[!is.na(coding_win$best)]

coding_win$best_model <- coding_win$braker_best
coding_win$best_model[is.na(coding_win$braker1_gene) & is.na(coding_win$braker2_gene) & !is.na(coding_win$stringtie_gene) ] <- 
  coding_win$stringtie_gene[is.na(coding_win$braker1_gene) & is.na(coding_win$braker2_gene) & !is.na(coding_win$stringtie_gene) ]

stringtie_best <- findBestHomolog(coding_win=coding_win,homolog_mode="stringtie")

coding_win<- left_join(coding_win,stringtie_best,by=c("best_model"="worst"))

coding_win$best_model[!is.na(coding_win$best.y)] <- coding_win$best.y[!is.na(coding_win$best.y)]

busco_convert <- makeBuscoConvert(coding_win=coding_win,genome_busco=genome_busco)
  
coding_win_busco <- left_join(coding_win,busco_convert,by=c("best_model"="value"))

coding_win_busco$final_gene <- coding_win_busco$best_model
coding_win_busco$final_gene[!is.na(coding_win_busco$dupName)] <- coding_win_busco$dupName[!is.na(coding_win_busco$dupName)]

#
coding_win_busco$source <- NA
coding_win_busco$source[coding_win_busco$best_model %in% na.omit(unique(coding_win_busco$braker1_gene))] <- "BRAKER1"
coding_win_busco$source[coding_win_busco$best_model %in% na.omit(unique(coding_win_busco$braker2_gene))] <- "BRAKER2"
coding_win_busco$source[coding_win_busco$best_model %in% na.omit(unique(coding_win_busco$stringtie_gene))] <- "stringtie"

noncoding <- addNonCoding2Stringtie(CPCfile='genome_files/annotation_files/FSJ2.stringtie_cpc.txt')
coding_win_nc <- left_join(coding_win_busco,noncoding)

save(coding_win_nc,file="genome_files/FSJ2.annotation_win.rdata")


####make new annotation .gtf file####                     

braker1_filt <- filterGTF(
  predictionType = "braker1",
inputFile='genome_files/annotation_files/FSJ2.BRAKER1.gtf.gz',
gene_list = finalGene_tally$best_model[finalGene_tally$source == "BRAKER1"])

braker1_filt$V9 <- gsub("file_1","b1",braker1_filt$V9)

braker2_filt <- filterGTF(
  predictionType = "braker2",
  inputFile='genome_files/annotation_files/FSJ2.BRAKER2.gtf.gz',
  gene_list = finalGene_tally$best_model[finalGene_tally$source == "BRAKER2"])

braker2_filt$V9 <- gsub("file_1","b2",braker2_filt$V9)


stringer_filt <- filterGTF(
  predictionType = "stringtie",
  inputFile= 'genome_files/annotation_files/FSJ2.stringtieMerged.gtf',
  gene_list = finalGene_tally$best_model[finalGene_tally$source == "stringtie"])

genes_filt <- 
  rbind.data.frame(braker2_filt,braker1_filt)

write.table(genes_filt, file = paste("genome_files/annotation_files/FSJ2.annotation.",today,".gtf",sep=""), 
            append = F, quote = FALSE, sep = "\t", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))
