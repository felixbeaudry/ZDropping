setwd('~/Documents/Github/genomeContent/')

options(scipen=999)

library(tidyverse)
library(ff)

library(data.table)
library(ggplot2)
library(cowplot)
library(FactoMineR) #PCA

library(qualV)
require(splines)

library(foreach)
library(MASS) #ditdistr

source("linkerScripts/linkedsel_functions.R")
est_rec_piece_func  <- function(x) {est_rec_piece(x)}

'%ni%' <- function(x,y)!('%in%'(x,y))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

assemblySizes <- function(sizesFile){  
 sizes <- sizesFile
  names(sizes) <- c("scaff","bp")

  sizes <- sizes[order(sizes$bp),]
  sizes$sizescum <-  cumsum(sizes$bp)
  
  N50 <-  unlist(sizes[sizes$sizescum - (sum(sizes$bp)/2) >= 0,][1,2])
  
  cat("Number of scaffolds is: ",length(sizes$bp),'\n')
  cat("genome size is: ",sum(sizes$bp),'bp (',sum(sizes$bp)/1000000000,'gb)\n')
  cat("genome N50 is: ",N50[1],"bp (",N50[1]/1000000,"mb)\n")
  
  print(ggplot(sizes,aes(x=bp/1000000)) + geom_histogram(bins=30))
  return(sizes[,-3])
}

windowedMap <- function(window_size,sizes){
  
  m_map <-
    foreach(chrom1=unique(sizes$scaff[sizes$bp > window_size]),.combine=rbind) %do% {
      win_tmp <-foreach(win1=seq(0,sizes$bp[sizes$scaff == chrom1],window_size),.combine=rbind) %do% {
        pos_tmp <- c(chrom1,win1)
        
      }
      win_tmp
    }
  
  df_map <- as.data.frame(m_map)
  df_map$V2 <- as.numeric(as.character(df_map$V2))
  df_map$mat_pos <- seq(1,length(df_map$V2 ),1)
  df_map$mb <- df_map$V2 / 1000000
  return(df_map)
}

homologyMap <- function(conversionFile,sizes){
  
  homology_conversion <- fread(conversionFile)
  names(homology_conversion) <- c("homolog","scaffold")
  
  size_homology <- left_join(sizes,homology_conversion,by=c("scaff"="scaffold"))
  
  size_homology$homolog_up <- factor(size_homology$homolog, levels = c(unique(na.omit(size_homology$homolog[order(-size_homology$bp)]))))
  size_homology$homolog_down <- factor(size_homology$homolog, levels = c(unique(na.omit(size_homology$homolog[order(size_homology$bp)]))));
  
  #stats on scaffolds with/without homologs
  cat("Homology search returned",length(unique(size_homology$scaff[!is.na(size_homology$homolog)])),"scaffolds with significant homology to",length(unique(size_homology$homolog)),"chromosomes\n" )
  
  cat("\n\nAssembly stats for scaffolds with homology\n")
  assemblySizes(size_homology[!is.na(size_homology$homolog),c(1,2)])
  cat("\n\nAssembly stats for scaffolds without homology\n")
  assemblySizes(size_homology[is.na(size_homology$homolog),c(1,2)])
  
  return(size_homology)
}

makeGeneticMap <- function(genmapFile,badMarkerDistance=30){
  genmap <- fread(genmapFile,fill=TRUE)
  
  genmap <- genmap[,c(1:6)]
  names(genmap) <- c('lg','raworder','marker','cm','V5','distance2Next')
  
  #remove faulty markers 
  genmap_filt <- genmap %>% filter(distance2Next < badMarkerDistance)
  genmap_filt <- genmap_filt %>% group_by(lg) %>% mutate(cm_adj= cm - min(cm))
  
  markerPos_plot<- 
    ggplot(genmap_filt ,aes(x=as.factor(lg),y=cm_adj)) +
    geom_line() + 
    geom_point(shape=3)+ 
    theme_bw(base_size=12) +
    labs(x="Linkage group", y="Genetic Distance (cM)") +
    theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    theme(panel.spacing = unit(0.1, "lines")) 
  
  
  #map length
  totalLength = 0
  for(chr in unique(genmap_filt$lg)){
    totalLength = max(genmap_filt$cm[genmap_filt$lg == chr]) + totalLength
  }
  
  genmap_tally <- as.data.frame(genmap_filt %>% group_by(lg) %>% tally())
  
  markerCount_plot <- ggplot(genmap_tally,aes(x=as.factor(lg),y=n)) + 
    geom_bar(stat="identity") + theme_bw(base_size = 12) + labs(x="Linkage group",y="Markers (SNPs)")
  
  print(plot_grid( markerCount_plot,  markerPos_plot,  labels = c(' ', ' '), ncol = 1, align = 'vh',axis='tbrl'))
  
  
  cat( 
    "\nNumber of linkage groups: ",length(unique(genmap_filt$lg)) ,
    "\nTotal Number of markers: ",length(unique(genmap_filt$marker)),
    "\nTotal map length: ",totalLength,
    "\n")
  
  return(genmap_filt)
  
}

makeMareyMap <- function(chipFile,genmap,size_homology){
  chip <- fread(chipFile)
  
  marey <- dplyr::left_join(chip,genmap,by=c("SNPname"="marker"))
  
  marey$mb <- as.numeric(marey$SNPpos) / 1000000
  
  marey$NewScaff[marey$NewScaff == "*"] <- "Unaligned"
  
  ##add chromosome numbers, from sequence homology
  marey <- left_join(marey,size_homology,by=c("NewScaff"="scaff"))
  
  marey$homolog[is.na(marey$homolog) ] <- "Unplaced"
  
  ##flip chromosome SNP mb positionss so they follow the marey map
  marey_zerod <-  marey[order(marey$homolog,marey$cm),]
  marey_zerod$cm_flip <- marey_zerod$cm
  
  for(chrom in unique(marey_zerod$homolog[marey_zerod$homolog != "Unplaced"])){
    marey_temp <- marey_zerod[marey_zerod$homolog == chrom,]
    if(length(unique(marey_temp$cm))>4){
      if(cor.test(marey_temp$cm,marey_temp$mb)[4] < 0){
        marey_zerod$cm_flip[marey_zerod$homolog == chrom] <- 
          max(marey_zerod$cm[marey_zerod$homolog == chrom],na.rm = T) - marey_zerod$cm[marey_zerod$homolog == chrom]
      }
    }
  }
  
  marey_zerod <- marey_zerod %>% filter(!is.na(cm_flip) & !is.na(mb) )
  
  marey_max <-  marey_zerod %>% group_by(homolog) %>% summarise(max_cm = max(cm_flip), max_mb=max(mb))
  marey_max$cmmb <- marey_max$max_cm/marey_max$max_mb
  marey_max<-  marey_max[complete.cases(marey_max),]
  marey_max$max_log_mb <- log(marey_max$max_mb)
  
  
  fullMap_plot <- 
    ggplot(marey_zerod %>% filter(homolog != "Unplaced") ,aes(x=mb,y=cm_flip)) +
    
    geom_point(shape=3,alpha=0.5)+ 
    
    facet_wrap(. ~ homolog, scales = 'free',  strip.position="bottom") +
    labs(x="", y=" ",fill=" ") +
    theme_bw(base_size=12) + 
    guides(fill="none")+
    theme(strip.background =element_rect(fill="white"),axis.title.x=element_blank()) + #xlim(120,125) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    theme(panel.spacing = unit(0.1, "cm"),
          panel.border = element_rect(color = "white", fill = NA, size = 0)) 
  
  placedMarkers_plot <- 
    ggplot(marey ,aes(x=as.factor(lg),y=cm)) +
    
    geom_line() + 
    geom_point(aes(shape=ifelse(is.na(homolog),"Unmapped","Mapped")))+ 
    scale_shape_manual(name="",values=c(Unmapped=19,Mapped=3))+
    
    theme_bw(base_size=12) +
    
    labs(x="LG", y="Genetic Distance (cM)") +
    theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    theme(panel.spacing = unit(0.1, "lines")) 
  
  cm2mb_plot <- 
    ggplot(marey_max,aes(y=log(cmmb),x=-max_log_mb)) +
    stat_smooth(method="lm",color="#FFB14E",formula= 'y ~ x') +
    
    geom_point(size=6) +
    geom_text(aes(label=homolog),color="white") + 
    theme_bw(base_size=16) +
    
    labs(x="-log Chromosome Length (Mb)",y="log Recombination Intensity (cM/Mb)") +
    theme(aspect.ratio = 1)
  
  cat("plot 1/3: Mapped Markers\n")
  suppressWarnings(print(placedMarkers_plot))
  
  cat("plot 2/3: cM to mb genome-wide\n")
  suppressWarnings(print(cm2mb_plot))
  
  cat("plot 3/3: all Marey Maps\n")
  suppressWarnings(print( fullMap_plot))
  
  
  cat(
    "Number of mapped markers: ", length(marey$SNPname[marey$homolog != "Unplaced"]),
    "\nNumber of scaffolded markers: ",length(marey$SNPname[marey$homolog != "Unplaced" & !is.na(marey$cm)]) 
  )
  
  return(marey_zerod)
}

get_rec_windows <- function(CMMB,printPlot=T,wSize = 1){
  names(CMMB) <- c("chr", "marker", "cm", "mb")
  CMMB <- CMMB[complete.cases(CMMB),]
  CMMB <- CMMB[order(CMMB$chr,CMMB$cm),]
  
  CMMB_reform <- cbind.data.frame(CMMB,"sp"="FSJ")
  CMMB_clean <- cleanup_maps(CMMB_reform)
  CMMB_mmRem <- remove_mismapped(CMMB_clean) # remove mismapped markers from all_maps_clean
  
  rec_piecewise <-  est_rec_piece_func(CMMB_mmRem[CMMB_mmRem$good.marker == TRUE, c("chr", "cm", "mb"),])
  
  #estimate rec rate at each SNP based on model
  CMMB$chrchr <- paste("chr",CMMB$chr,sep="")
  
  
  window.size = wSize 
  
  loop = 0
  for(z in unique(CMMB$chrchr)){
    
    mb <- sort(CMMB$mb[CMMB$chrchr == z])
    
    Y.piece.data = data.frame(mb=mb-(window.size/2000))
    Y.piece.data = rbind(Y.piece.data, tail(mb, n=1)+(window.size/2000))
    
    #predict
    piece.fun=rec_piecewise[[z]][["func"]]
    Y.piece.raw = predict(piece.fun, newdata=Y.piece.data)
    
    #make table
    Y.piece.raw[Y.piece.raw < 0] = 0 #no negative cM
    
    if(loop == 0){
      loc_scaf_cM <-  cbind.data.frame("chrchr"=z,"mb"=mb, "cm" = Y.piece.raw[-c(1)])
      
    }else{
      loc_scaf_cM <- rbind.data.frame(loc_scaf_cM,
                                      cbind.data.frame( "chrchr"=z, "mb"=mb, "cm" = Y.piece.raw[-c(1)])
      )
    }
    loop =+ 1
  }
  
  fitVraw <- 
    rbind.data.frame(
      cbind.data.frame(CMMB[,c(5,4,3)],fit="raw"),
      cbind.data.frame(loc_scaf_cM,fit="fit")
    )
  
  if(printPlot==T){
    print(ggplot(fitVraw,aes(x=mb,y=cm,color=fit)) + geom_point(size=0.1)  +  facet_wrap(. ~ chrchr,scales = "free") + theme_bw())
  }
  
  #calculate rate rec in windows
  windowedrec <- cbind.data.frame("LG"=NA,"winNum"=NA,"mb_s"=NA,"mb_e"=NA,"cmmb"=NA)
  
  for(z in unique(loc_scaf_cM$chrchr)){
    loc_scaf_cM_loop <- loc_scaf_cM[loc_scaf_cM$chrchr == z,]
    mb <- seq(from=0,to=max(loc_scaf_cM_loop$mb),0.001)
    
    Y.piece.data = data.frame(mb=mb-(window.size/2000))
    Y.piece.data = rbind(Y.piece.data, tail(mb, n=1)+(window.size/2000))
    
    #predict
    piece.fun=rec_piecewise[[z]][["func"]]
    Y.piece.raw = predict(piece.fun, newdata=Y.piece.data)
    
    #make table
    Y.piece.raw[Y.piece.raw < 0] = 0 #no negative cM
    loc_scaf_cM_loop <-   cbind.data.frame(  mb, "cM" = Y.piece.raw[-c(1)])
    
    win = 0
    wl = 0
    
    while(wl < max(loc_scaf_cM_loop$mb)){
      wh = wl + wSize
      
      loc_scaf_cM_loop_w <- loc_scaf_cM_loop[loc_scaf_cM_loop$mb > wl & loc_scaf_cM_loop$mb <= wh,]
      
      if(all(is.na(loc_scaf_cM_loop_w$mb)) || all(is.na(loc_scaf_cM_loop_w$cM))){
        localRec <- NA
        
      } else{
        localRec <- lm(loc_scaf_cM_loop_w$cM ~ loc_scaf_cM_loop_w$mb  ) 
        localRec <- localRec$coefficients[2]
        if(localRec < 0 || is.na(localRec)){
          localRec <- NA
        }
      }
      
      windowedrec <- rbind.data.frame(windowedrec, 
                                      cbind.data.frame("LG"=z,"winNum"=win,"mb_s"=wl,"mb_e"=wh,"cmmb"=localRec)
      )
      win = win + 1
      wl = wl + wSize
      
    }
    
  }
  
  windowedrec <- separate(windowedrec, LG, c("ch","LG"), 
                          sep = "r", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
  
  windowedrec <- windowedrec[complete.cases(windowedrec), c(2,5,6)]
  
  return(windowedrec)
}

ldhatFit <- function(ldhatFile,size_homology){
  
  load(ldhatFile)
  LD.list.ordered_sub$mb <- LD.list.ordered_sub$POS / 1000000
  LD.list.ordered_sub$Cum_rho_adj <- LD.list.ordered_sub$Cum_rho/5000 #approx Ne?
  
  LD.list.ordered_sub_joined <-left_join(LD.list.ordered_sub,size_homology,by=c("CHROM"="scaff"))

  return(LD.list.ordered_sub_joined)
}

compareRHO2MB <- function(windowedreC_join,size_homology){
  windowedreC_comp <- windowedreC_join[complete.cases(windowedreC_join),]
  
  #filter
  windowedreC_filt <- windowedreC_comp %>% filter(rho > 2 | link < 6 | xor(rho < 2 , link > 6))
  
  #calculate correlations for each window
  rec_corr <- foreach(chrom=unique(windowedreC_join$LG),.combine=rbind) %do% {
    windowedreC_tmp <- windowedreC_filt %>% filter(LG == chrom)
    rec_corr_tmp <- c(chrom,cor(windowedreC_tmp$rho,windowedreC_tmp$link,method = "spearman"))
    rec_corr_tmp
  }
  rec_corr <- as.data.frame(rec_corr)
  names(rec_corr) <- c("chrom","corr")
  rec_corr$corr <- as.numeric(as.character(rec_corr$corr))
  
  windowedreC_corr <- left_join(windowedreC_filt,rec_corr,by=c("LG"="chrom")) %>% left_join(size_homology,by=c("LG"="homolog")) 

  print(
    ggplot(windowedreC_corr,aes(x=link,y=rho)) + 
      geom_point(alpha=0.75,size=0.75) + 
      stat_smooth(method="lm",formula = 'y~x',se=FALSE,aes(color=corr))+  
      facet_wrap(. ~homolog_up,scales="free") + theme_bw() +
      labs(x="Linkage (cM/mb)",y="LD (rho)",color="Corr. (Sp.Rank)") + 
      scale_color_gradient(
        low = "white",
        high = "black",
        na.value = "gray"
      ) 
  )
  return(windowedreC_corr)
}

makeSexMarey <- function(sexlinkagefile,marey){
  
  sex_map_raw <-  fread(sexlinkagefile)[,-10]
  names(sex_map_raw) <- c("lg","markerNum","marker","male","female","male_r","male_r_lnl","female_r","female_r_lnl")
  
  #remove unlikely markers
  sex_map <- sex_map_raw %>% filter(marker %in% marey$SNPname) 
  sex_map <- sex_map %>% group_by(lg) %>% mutate(male_adj= male - min(male),female_adj= female - min(female))
  
  sex_max <- sex_map %>% group_by(lg) %>% summarise(maxmale = max(male),maxfem = max(female))
  
  sex_dist_plot <- 
    ggplot(sex_max,aes(x=-log(maxmale),y=-log(maxfem))) +
    geom_abline(intercept = 0, slope = 1,linetype="dashed",color="grey")+
    geom_smooth(method="lm",color="#FFB14E") +
    geom_point(size=6)+
    geom_text(aes(label=lg),color="white") +
    theme_bw(base_size=16) + 
    labs(x="-log Male Linkage Map length (cM)",y="-log Female Linkage Map length (cM)") +
    theme(aspect.ratio = 1)
  
  sex_map_melt <- melt(sex_map %>% dplyr::select(lg,marker,male_adj,female_adj) ,id.vars = c("marker","lg"),verbose=FALSE) #[,c(1,3,10:11)]
  
  sex_map_melt$variable_name[sex_map_melt$variable == "male_adj"] <- "male"
  sex_map_melt$variable_name[sex_map_melt$variable == "female_adj"] <- "female"
  
  sex_map_plot <-
    ggplot(sex_map_melt ,aes(x=lg,y=value,color=variable_name)) +
    
    geom_line(position=position_dodge(width=0.5)) + 
    geom_point(position=position_dodge(width=0.5),shape=3)+ 
    
    theme_bw(base_size=24) +
    
    labs(x="LG", y="Genetic Distance (cM)",color="") +
    theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    theme(panel.spacing = unit(0.1, "lines"))  + theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    ) + scale_color_manual(values = c( "#00ADEF", "#00A550"))
  
  #print out longer map for each lg
  totalLength_m = 0
  totalLength_f = 0
  for(chr in unique(sex_map_melt$lg)[-1]){
    totalLength_m = max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable_name == "male"],na.rm=T) + totalLength_m
    totalLength_f = max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable_name == "female"],na.rm=T) + totalLength_f
    if(max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable_name == "male"],na.rm=T) > max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable_name == "female"],na.rm=T)){cat("LG",chr,": male map longer\n")}
    if(max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable_name == "male"],na.rm=T) < max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable_name == "female"],na.rm=T)){cat("LG",chr,": female map longer\n")}
    
  }
  
  ##with mb position
  sex_marey <- left_join(sex_map_melt,marey %>% dplyr::select(SNPname,mb,cm_adj,homolog,homolog_up),by=c("marker"="SNPname")) #[,c(1,13,12)]
  
  ##flip chromosome SNP mb positionss so they follow the marey map
  sex_marey_zerod <-  sex_marey[order(sex_marey$lg,sex_marey$variable,sex_marey$value),]
  sex_marey_zerod$cm_flip <- sex_marey_zerod$value
  
  for(chrom in unique(sex_marey_zerod$lg)){
    # cat(chrom)
    marey_temp <- sex_marey_zerod[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name == "male",] #chose male so flip can work on Z too
    
    if(cor.test(marey_temp$value,marey_temp$mb)[4] < 0){
      sex_marey_zerod$cm_flip[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "male"] <- 
        max(sex_marey_zerod$value[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "male"]) - sex_marey_zerod$value[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "male"]
      
      sex_marey_zerod$cm_flip[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "female"] <- 
        max(sex_marey_zerod$value[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "female"]) - sex_marey_zerod$value[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "female"]
    }
  }
  
  sex_marey_zerod <- sex_marey_zerod[complete.cases(sex_marey_zerod),]
  
  fullSexMap_plot <- 
    ggplot(sex_marey_zerod  ,aes(x=mb,y=cm_flip,color=variable_name)) +
    
    geom_point(shape=3,alpha=0.5)+ 
    
    facet_wrap(. ~ homolog_up, scales = 'free',  strip.position="bottom") +
    labs(x="", y=" ",fill=" ") +
    theme_bw(base_size=12) + 
    guides(fill="none")+
    theme(strip.background =element_rect(fill="white"),axis.title.x=element_blank()) + #xlim(120,125) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    theme(panel.spacing = unit(0.1, "cm"),
          panel.border = element_rect(color = "white", fill = NA, size = 0))  + 
      scale_color_manual(values = c( "#00ADEF", "#00A550"))
  
  
  cat("Plot 1/3: Male to Female Recombination Rate\n")
  print(sex_dist_plot)
  cat("Plot 2/3: Male and Female Genetic Distances\n")
  print(sex_map_plot)
  cat("Plot 3/3: Male and Female Marey Maps\n")
  print(fullSexMap_plot)
  return(sex_marey_zerod)
}

makeLinkageSDI <- function(rec_windows_male,rec_windows_female,marey){
  names(rec_windows_male)[3] <- "male"
  names(rec_windows_female)[3] <- "female"
  
  marey$lg <-  as.character( marey$lg)
  
  rec_windows_sex <- left_join(rec_windows_male,rec_windows_female) %>% left_join(marey %>% dplyr::select(lg,homolog),by=c("LG"="lg"))
  
  rec_windows_sex$male <- as.numeric(rec_windows_sex$male)
  rec_windows_sex$female <- as.numeric(rec_windows_sex$female)
  
  rec_windows_sex$SDI <- ifelse(rec_windows_sex$male >rec_windows_sex$female,-1*((rec_windows_sex$male/rec_windows_sex$female) - 1), 1*((rec_windows_sex$female/rec_windows_sex$male) - 1))

  return(rec_windows_sex)
} 

calcGapsNBases <- function(bybaseFile,homology_conversion,gap_window){
  
  bybase <- fread(bybaseFile)
  bybase$scaffold <- gsub(pattern = ">","",bybase$scaffold)
  
  #remove last window on each scaffold (they are less than full length)
  bybase_filt <- bybase %>% group_by(scaffold) %>% filter(bp != max(bp)) 
  
  bybase_filt <- left_join(bybase_filt,homology_conversion,by=c("scaffold"="scaff"))
  
  
  bybase_filt$N_frac <- bybase_filt$N / gap_window
  
  bybase_filt$GC <- (bybase_filt$G + bybase_filt$C)/ gap_window
  bybase_filt$AT <- (bybase_filt$A + bybase_filt$T)/ gap_window
  bybase_filt$mb <- bybase_filt$bp.x / 1000000
  
  
  #names(bybase) #"chrom"  "mb"     "N_frac"      "GC"   "AT"
  print(
    ggplot(bybase_filt %>% filter(!is.na(homolog_up)),aes(x=homolog_up,y=GC)) + 
      geom_boxplot() +
      theme_bw(base_size=15) + guides(fill="none") + 
      theme(strip.background =element_rect(fill="white"),axis.title.x=element_blank()) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_fill_manual(values = c(
        "cornflowerblue",
        "indianred"))
  )
  

  
 # bybase_melt <- reshape2::melt(bybase_filt %>% dplyr::select(scaffold,homolog,mb,N_frac,  GC , AT),id.vars = c("homolog","mb","scaffold"),verbose=FALSE)
 # bybase_melt$type <- factor(bybase_melt$variable, levels=c( "GC", "AT", "N_frac"))
  
  return(bybase_filt %>% dplyr::select(scaffold,homolog,mb,N_frac,  GC , AT))
}

makeCpG <- function(cpgFile,homology_conversion){
  
  
  cpg_tally <- fread(cpgFile) %>% group_by(chrom) %>% tally() %>% left_join(homology_conversion,by=c("chrom"="scaff"))
  
  cpg_tally$cpgMb <- cpg_tally$n / (cpg_tally$bp/1000000)
  
  names(cpg_tally)[2] <- "CpG"
  
  print( 
    ggplot(cpg_tally,aes(x=homolog_up,y=cpgMb)) + 
      geom_boxplot() + 
      theme_bw(base_size=15) + guides(color="none") + 
      theme(strip.background =element_rect(fill="white"),axis.title.x=element_blank()) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      scale_color_manual(values = c(
        "cornflowerblue",
        "indianred")) + labs(y="CpG Islands/mb")
  )
  return(cpg_tally)
}

makeRepeatWindows <- function(genomeName,repeatFile=NA, window_size_TE ,df_map, homology_conversion){
  if(is.na(repeatFile)){
    repeatFile=paste0('genome_files/',genomeName,'.repeats.txt.gz')
  }
  All_repeats <- fread(repeatFile,header=FALSE, sep=" ",fill = TRUE,skip=3) 
  names(All_repeats) <- c("SWscore","pdiverged","pdeleted","pinserted","scaff","LG_start","LG_end","left","dir","repeatname","family","rep_start","rep_end","end_left","ID","star")
  
  All_repeats <- separate(All_repeats, family,  c("family","subfam"), 
                          sep = "/", remove = TRUE, convert = FALSE, extra = "merge", fill = "right")
  
  #filter and re-categorize
  All_repeats <- All_repeats[All_repeats$family != "ARTEFACT",] 
  All_repeats$family[All_repeats$family %in% c("Unknown","Unspecified","Other")] <- "Unspecified"
  All_repeats$family[All_repeats$family %in% c("LTR?","LTR")] <- "LTR"
  All_repeats$family[All_repeats$family %in% c("SINE","SINE2","SINE?")] <- "Non-LTR"
  All_repeats$family[All_repeats$family %in% c("LINE","LINE?")] <- "LINE"
  All_repeats$family[All_repeats$family %in% c("DNA","Tc1-Mariner","RC","RC?","DNA?")] <- "DNA_transposon"
  All_repeats$family[All_repeats$family %in% c("tRNA","rRNA","snRNA","scRNA")] <- "house_RNAs"
  All_repeats$family[All_repeats$family %in% c("Low_complexity")] <- "Simple_repeat"
  
  All_repeats_mut <-
    All_repeats %>%
    mutate(LG_lower_start=pmin(LG_start,LG_end), featuresize=abs(LG_start-LG_end)) %>%
    group_by(scaff) %>%
    mutate( LG_lower_start_adj=LG_lower_start-min(LG_lower_start))
  
  cat("Total repeat length:",sum(All_repeats_mut$featuresize)/1000000 ,"mb\n")
  
  cat("repeats account for",(sum(All_repeats_mut$featuresize) / sum(homology_conversion$bp))*100,"% of the assembly\n")
  
  All_repeats_lengths <- as.data.frame(All_repeats_mut %>% group_by(family) %>% summarise(length=sum(featuresize),n=n()))
  All_repeats_lengths$frac <-   All_repeats_lengths$length / sum(homology_conversion$bp)
  
  cat("The most common element is...\n")
  print(All_repeats_lengths[order(-All_repeats_lengths$n),][c(1:2),])
  
  cat("The element making up the largest fraction of the genome is...\n")
  print(All_repeats_lengths[order(-All_repeats_lengths$frac),][c(1:2),])
  
  All_repeats_win <-  
    All_repeats_mut %>% group_by(scaff,family)  %>% 
    mutate(position_window=LG_lower_start_adj%/%as.numeric(window_size_TE)) %>% 
    group_by(scaff,position_window,family) %>%  
    add_tally() %>% select_if(., is.numeric) %>% 
    summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) %>% 
    mutate(proportion_of_window=featuresize_sum/as.numeric(window_size_TE)) %>% 
    mutate(window_start=(position_window*(as.numeric(window_size_TE)+1)) )
  
  All_repeats_win$mb <- round(as.numeric(All_repeats_win$window_start) / 1000000)
  
  ##add in empty windows for each family: make data.frame with full length of windows for each chrom/scaff 
  
  df_map_TE <- 
    foreach(TEfam=unique(All_repeats_win$family),.combine=rbind.data.frame) %do% {
      cbind.data.frame(df_map,"TEfam"=TEfam )
    }
  
  All_repeats_zeros <- left_join(df_map_TE,All_repeats_win,by=c("V1"="scaff","mb"="mb","TEfam"="family"))
  
  All_repeats_zeros$proportion_of_window[is.na(All_repeats_zeros$proportion_of_window)] <- 0
  
  All_repeats_zeros <- left_join(All_repeats_zeros,homology_conversion,by=c("V1"="scaff"))
  
  #pdf(paste0("images/",genomeName,".repeats.pdf"),width=4.5,height=3)
  All_repeats_zeros_filt <- All_repeats_zeros %>% filter(!is.na(homolog_up))
  
  print( 
    
    ggplot( All_repeats_zeros_filt,aes(x=mb,y=proportion_of_window,fill=TEfam)) + 
      
      geom_area(position="stack") + 
      
      facet_wrap(. ~ homolog_up, scales = 'free',strip.position = "bottom") +
      labs(x="", y=" ",fill=" ") +
      theme_bw(base_size=8) + 
      theme(strip.background =element_rect(fill="white"),axis.title.x=element_blank()) + #ylim(0,1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      theme(panel.spacing = unit(-0.1, "cm"),
            panel.border = element_rect(color = "white", fill = NA, size = 0))  +
      scale_fill_manual(values=c("#ffd700",
                                 "#ffb14e",
                                 "#fa8775",
                                 "#ea5f94",
                                 "#cd34b5",
                                 "#9d02d7",
                                 "#0000ff","grey"))
  )
  # dev.off()
  
  return(All_repeats_zeros)
}

hic_summarizer <- function(genomeName,window_size,homology_conversion){
  hic_raw <- fread(paste0('genome_files/',genomeName,'.hic.',window_size/1000000,'mb.raw.txt'))
  
  hic_m_map <- left_join(hic_raw,homology_conversion,by=c("V1"="scaff")) %>% left_join(homology_conversion,by=c("V2"="scaff"))
  
  
  hic_inter <- hic_m_map %>% group_by(homolog_up.x,homolog_down.y) %>% 
    summarise(mean_connect = mean(V5))
  
  hic_filt <- hic_inter %>% 
    filter(homolog_up.x != homolog_down.y & !is.na(homolog_up.x) & !is.na(homolog_down.y))
  
  
  print(
    ggplot(hic_filt, aes(x=homolog_up.x, y=homolog_down.y, fill= log(mean_connect))) + 
      geom_tile(color="lightgray") + #theme_bw(base_size=18) + 
      scale_fill_gradient2(
        low = "cornflowerblue",
        mid = "white",
        high = "indianred",
        midpoint = mean(log(hic_filt$mean_connect)),
        na.value = "darkgray"
      ) + labs(y="",fill="Hi-C links") + #guides(fill=FALSE) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank())   +
      theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5))
  )
  
  names(hic_filt) <- c("chr_1","chr_2","meanHiC")
  
  return(hic_filt)
}

coding_windower <- function(annotation_winData, size_homology,window_size,genomeName,busco_N = 8338){
  
  genome_busco <- getGenomeBUSCO(genomeBUSCOFile=paste0('genome_files/annotation_files/',genomeName,'.genome_busco.tsv'))
  
  load(annotation_winData)
  
  coding_win_nc$label[is.na(coding_win_nc$final_gene)] <- NA
  coding_win_nc$stringtie_bp[coding_win_nc$stringtie_bp > 1000] <- 1000
  
  unique(coding_win_nc$source)
  
  coding_win_nc$final_bp <- ifelse(coding_win_nc$source == "stringtie",
                                   coding_win_nc$stringtie_bp,
                                   ifelse(coding_win_nc$source == "BRAKER1",
                                          coding_win_nc$braker1_bp,
                                          ifelse(coding_win_nc$source == "BRAKER21",
                                                 coding_win_nc$braker2_bp,
                                                 0)
                                   )
  )
  
  cat("number of unique annotations:",length(unique(coding_win_nc$final_gene)) ,"genes\n")
  cat("number of unique protein-coding genes:",length(unique(coding_win_nc$final_gene[coding_win_nc$label == "coding"])) ,"genes\n")
  
  cat("BUSCO score of annotation:", (length(unique(c(coding_win_nc$braker1_busco,coding_win_nc$braker2_busco,coding_win_nc$stringtie_busco))) /busco_N)*100,"%\n")
 
  coding_win_nc$win_r <- round(coding_win_nc$win_s,digits = -log(window_size,10))
  
  coding_oldChrom <- coding_win_nc %>% group_by(chrom,win_r,label) %>% summarise(sumBP=sum(final_bp,na.rm = T))
  
  coding <- left_join(coding_oldChrom,size_homology,by=c("chrom"="scaff"))
  
  coding$mb <- coding$win_r/1000000
  coding$fracCoding <- coding$sumBP/window_size
  coding <-  coding %>% filter(!is.na(label))
  
  print( 
    
    ggplot( coding,aes(x=mb,y=fracCoding,fill=label)) + 
      
      geom_area(position="stack") + 
      
      facet_wrap(. ~ homolog_up, scales = 'free_x',strip.position = "bottom") +
      labs(x="", y=" ",fill=" ") +
      theme_bw(base_size=8) + 
      theme(strip.background =element_rect(fill="white"),axis.title.x=element_blank()) + #ylim(0,1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      theme(panel.spacing = unit(-0.1, "cm"),
            panel.border = element_rect(color = "white", fill = NA, size = 0))  +
      scale_fill_manual(values=c("#ffd700",
                                 "#ffb14e",
                                 "#fa8775",
                                 "#ea5f94",
                                 "#cd34b5",
                                 "#9d02d7",
                                 "#0000ff","grey"))
  )
  
  
  return(coding)
  
}

'%ni%' <- function(x,y)!('%in%'(x,y))

getGenomeBUSCO  <- function(genomeBUSCOFile,aves_busco_N = 8338){
  genome_busco <- fread(genomeBUSCOFile,sep="\t",header=TRUE,select = c(1:10),fill=TRUE)
  names(genome_busco)[2] <- "genome"
  cat("Genome-assembly BUSCO completeness: ", 100*round(length(unique(genome_busco$Busco_id[genome_busco$genome != "Missing"]))/aves_busco_N,3),"%\n")
  # buscos_raw <- unique(genome_busco$Busco_id)
  for(id in unique(genome_busco$Busco_id)){
    genome_busco$dupCount[genome_busco$Busco_id == id] <- c(1:length(genome_busco$Busco_id[genome_busco$Busco_id == id]))
  }
  genome_busco$dupName <- paste(genome_busco$Busco_id,genome_busco$dupCount,sep="d")
  return(genome_busco)
  
  
}

reformatGTF <- function(inputFile,predictionType,buscoFile=NA,aves_busco_N = 8338){
  cat("## Reformatting ",predictionType,"##\n" )
  
  #autoload BUSCO file if BUSCO file name is not passed to function
  if(is.na(buscoFile)){
    buscoFile <- paste(predictionType,'.busco.tsv',sep="")
    cat(paste("**assuming busco file is called ",predictionType,'.busco.tsv**',sep=""))
    
  }
  #load BUSCO
  busco <- fread(buscoFile,sep="\t",header=TRUE,fill=TRUE)
  
  #split predictino type into stringtie and BRAKER1
  if(predictionType == "stringtie"){
    annotation <- fread(inputFile,sep = "\t", skip = 2)
    
    ann_filt <- annotation %>% filter(V3 == "exon") %>% arrange(V1, V4)
    split_gtf <-
      separate(ann_filt, V9, c("gene","gene_id","transcript","transcript_id","exon","exon_id","finalColumns"), 
               sep = "\"", remove = TRUE, convert = FALSE, extra = "merge", fill = "right")[,c(1:8,10,12)]
    busco <- separate(busco, Sequence, c("transcript_id","pos"), 
                      sep = "_stringtie:", remove = TRUE, convert = FALSE, extra = "merge", fill = "right")
    
    
  }
  
  if(predictionType %in% c("braker1","braker2")){
    annotation <- fread(inputFile,sep = "\t")
    
    ann_filt <- annotation %>% filter(V3 == "CDS") %>% arrange(V1, V4)
    split_gtf <-
      separate(ann_filt, V9, c("transcript","transcript_id","gene","gene_id","finalColumns"), 
               sep = "\"", remove = TRUE, convert = FALSE, extra = "merge", fill = "right") [,c(1:8,12,10)]
    
    busco <- separate(busco, Sequence, c("transcript_id","pos"),  sep = ":", remove = TRUE, convert = FALSE, extra = "merge", fill = "right")
    busco$transcript_id <- paste(predictionType,"_",busco$transcript_id,sep="")
    busco <- busco[,c(1,2,3,4,8,9,10,11)]
    
    split_gtf$gene_id <- gsub('file_1_file_1', predictionType, split_gtf$gene_id)
    split_gtf$transcript_id <- gsub('file_1_file_1', predictionType, split_gtf$transcript_id)
    
    
  }
  suppressMessages( exons <- left_join(split_gtf,busco) )
  
  names(exons) <- c('scaffold','origin','annotation','start_pos','end_pos',
                    'annotation_score','direction','readingFrame','gene_id','transcript_id',
                    'busco_id', 'busco_status','busco_pos','busco_score',
                    'busco_length','busco_gene_link','busco_gene_description')
  
  exons$transcript_id <- gsub('MSTRG.', 'STRG', exons$transcript_id)
  exons$gene_id       <- gsub('MSTRG.', 'STRG', exons$gene_id)
  
  cat("Unique Annotated Genes: ",length(unique(exons$gene_id)),"\n" )#unique genes
  
  cat("Annotation BUSCO completeness: ", round(length(unique(exons$busco_id)) / aves_busco_N,3)*100,"%","\n")
  exon_tally <- exons %>% group_by(transcript_id) %>% tally()
  cat("Unique Transcripts: ",sum(exon_tally$n) ,"\n")
  
  exons_transcript <- exons %>% group_by(gene_id,transcript_id) %>% tally()
  exons_gene <- exons_transcript  %>% group_by(gene_id) %>% tally()
  
  cat("Average Alternate Transcripts/gene: ",mean(exons_gene$n) ,"\n")
  cat("Max Alternate Transcripts/gene: ", max(exons_gene$n) ,"\n")
  
  exons$origin <- predictionType
  return(exons)
}

addNonCoding2Stringtie <- function(CPCfile){
  
  noncod_string <- fread(CPCfile)
  names(noncod_string)[1] <- "ID"
  noncod_string$ID <- gsub('MSTRG.', 'STRG',noncod_string$ID)
  noncod_string <- separate(noncod_string, ID, c("stringtie_gene","stringtie_t"),  sep = "\\.", remove = FALSE, convert = FALSE, extra = "merge", fill = "right")
  
  noncod_string_p <- noncod_string %>% group_by(stringtie_gene,label) %>% tally()
  noncod_string_p <- noncod_string_p[order(noncod_string_p$label, decreasing=FALSE),]
  noncod_string_p_clean <- noncod_string_p[!duplicated(noncod_string_p$stringtie_gene),]
  return(noncod_string_p_clean)
}

returnLongestTranscript <- function(genes_byGene,win,coding_window_size){
  longestTranscript <- 
    suppressMessages(
      genes_byGene %>% group_by(gene_id,transcript_id) %>% 
        summarise(codingLength = sum(gene_length),winStart=min(start_pos),winEnd=max(end_pos)) %>% 
        top_n(n=1,wt = codingLength)
    )
  longestTranscript <- longestTranscript[1,]
  
  if(longestTranscript$winStart < win){
    longestTranscript$codingLength <- longestTranscript$codingLength - (longestTranscript$winStart - win)
    
  }
  if(longestTranscript$winEnd > (win+coding_window_size)){
    longestTranscript$codingLength <- longestTranscript$codingLength - (longestTranscript$winEnd - (win+coding_window_size))
    
  }
  return( longestTranscript[,c(1,3)])
}

makeBuscoConvert <- function(coding_win,genome_busco,aves_busco_N = 8338){
  coding_win$win_sr <- round(coding_win$win_s,windowSizelog)
  genome_busco$win_sr <- round(genome_busco$Gene_Start,windowSizelog)
  
  braker1_buscoConvert <- coding_win %>% 
    filter( !is.na(braker1_busco) & braker1_gene %ni% c("bad_exon","overlapping_exons") ) %>% 
    dplyr::select(chrom,braker1_gene,braker1_busco,win_sr) %>% distinct()
  
  braker2_buscoConvert <- coding_win %>% 
    filter( !is.na(braker2_busco) & braker2_gene %ni% c("bad_exon","overlapping_exons") ) %>% 
    dplyr::select(chrom,braker2_gene,braker2_busco,win_sr) %>% distinct()
  
  stringtie_buscoConvert <- coding_win %>% 
    filter( !is.na(stringtie_busco) & stringtie_gene %ni% c("bad_exon","overlapping_exons") ) %>% 
    dplyr::select(chrom,stringtie_gene,stringtie_busco,win_sr) %>% distinct()
  
  all_busco <- 
    left_join(genome_busco,braker1_buscoConvert,by=c("Busco_id"="braker1_busco","Sequence"="chrom","win_sr"="win_sr")) %>%
    left_join(braker2_buscoConvert,by=c("Busco_id"="braker2_busco","Sequence"="chrom","win_sr"="win_sr")) %>%
    left_join(stringtie_buscoConvert,by=c("Busco_id"="stringtie_busco","Sequence"="chrom","win_sr"="win_sr"))
  
  cat("BUSCO completeness of combined annotation: ",
      (length(unique( all_busco$Busco_id[!is.na(all_busco$braker1_gene) |!is.na(all_busco$braker2_gene) |!is.na(all_busco$stringtie_gene)] ))/aves_busco_N)*100,
      "%")
  
  cat("\nBUSCO completeness compared to genome: ", 
      100*round(length(unique( all_busco$Busco_id[!is.na(all_busco$braker1_gene) |!is.na(all_busco$braker2_gene) |!is.na(all_busco$stringtie_gene)] ))/length(unique(genome_busco$Busco_id[genome_busco$genome != "Missing"])),3),"%\n")
  
  busco_convert <- all_busco %>% dplyr::select(dupName,braker1_gene,braker2_gene,stringtie_gene) %>% melt(id.vars="dupName") %>% filter(!is.na(value))
  
  return(busco_convert)
}

findBestHomolog <- function(coding_win,homolog_mode){
  if(homolog_mode == "braker"){
    Bests <- unique(coding_win$braker_best)
    best_tally <- coding_win %>% group_by(braker_best) %>% tally()
    
    homolog_inference <- coding_win %>% 
      filter(!is.na(braker1_gene) & !is.na(braker2_gene) & braker1_gene %ni% c("bad_exon","overlap") & braker2_gene %ni% c("bad_exon","overlap")) %>% 
      group_by(braker1_gene,braker2_gene) %>% tally() %>%
      arrange(braker1_gene,braker2_gene, -n) %>%
      filter(duplicated(braker1_gene) == FALSE)  %>% ungroup() %>%
      filter(duplicated(braker2_gene) == FALSE)
    
    homolog_best <- homolog_inference %>% dplyr::select("braker1_gene","braker2_gene") %>%
      left_join(best_tally,by=c("braker1_gene"="braker_best")) %>%
      left_join(best_tally,by=c("braker2_gene"="braker_best"))
  }
  if(homolog_mode == "stringtie"){
    Bests <- unique(coding_win$best_model)
    best_tally <- coding_win %>% group_by(best_model) %>% tally()
    
    homolog_inference <- coding_win %>% 
      filter(!is.na(braker_best) & !is.na(stringtie_gene) & stringtie_gene %ni% c("overlap")) %>% 
      group_by(braker_best,stringtie_gene) %>% tally() %>%
      arrange(braker_best,stringtie_gene, -n) %>%
      filter(duplicated(braker_best) == FALSE)  %>% ungroup() %>%
      filter(duplicated(stringtie_gene) == FALSE)
    
    homolog_best <- homolog_inference %>% dplyr::select("braker_best","stringtie_gene") %>%
      left_join(best_tally,by=c("braker_best"="best_model")) %>%
      left_join(best_tally,by=c("stringtie_gene"="best_model"))
    
    homolog_best$n.y <- 0
    
  }
  names(homolog_best)[c(1:2)] <- c("first","second")
  homolog_best$n.x[is.na(homolog_best$n.x)] <- 0
  homolog_best$n.y[is.na(homolog_best$n.y)] <- 0
  
  homolog_best$best <- ifelse(homolog_best$n.x >= homolog_best$n.y,
                              homolog_best$first,
                              homolog_best$second)
  homolog_best$worst <- ifelse(homolog_best$n.x >= homolog_best$n.y,
                               homolog_best$second,
                               homolog_best$first)
  
  best <- homolog_best %>% filter(worst %in% Bests) %>% ungroup() %>% dplyr::select(best,worst)
  return(best)
}

filterGTF <- function(inputFile,gene_list,predictionType){
  if(predictionType == "stringtie"){
    gtf_raw <- fread(inputFile,sep = "\t", skip = 2)
    gtf_sep <-
      separate(gtf_raw, V9, c("gene","gene_id","transcript","transcript_id","exon","exon_id","finalColumns"), 
               sep = "\"", remove = TRUE, convert = FALSE, extra = "merge", fill = "right")
    
    gtf_sep$row <- seq(1,length(gtf_raw$V1),by=1)
    
    gene_list_sub <- gsub('STRG','MSTRG.', gene_list)
    
    
    gtf_rows <- 
      sort(
        gtf_sep$row[gtf_sep$gene_id %in% gene_list_sub])
    
  }
  if(predictionType %in% c("braker1","braker2")){
    gtf_raw <- fread(inputFile)
    gtf_sep <-
      separate(gtf_raw, V9, c("transcript","transcript_id","gene","gene_id","end"), 
               sep = "\"", remove = TRUE, convert = FALSE, extra = "merge", fill = "right")
    
    gtf_sep$row <- seq(1,length(gtf_raw$V1),by=1)
    
    gtf_sep_g  <- gtf_sep %>% filter(V3 == "gene")
    gtf_sep_t <- gtf_sep %>% filter(V3 == "transcript")
    gtf_sep_t <- separate(gtf_sep_t, transcript, c("gene_id","transcript"),  sep = "\\.", remove = FALSE, convert = FALSE, extra = "merge", fill = "right")
    gtf_sep_o <- gtf_sep %>% filter(!is.na(gene_id) & V3 %ni% c("gene","transcript"))
    
    gene_list_sub <- gsub(predictionType,'file_1_file_1', gene_list)
    
    gtf_rows <- 
      sort(c(
        gtf_sep_g$row[gtf_sep_g$transcript %in% gene_list_sub],
        gtf_sep_o$row[gtf_sep_o$gene_id %in% gene_list_sub],
        gtf_sep_t$row[gtf_sep_t$gene_id %in% gene_list_sub]
      ))
  }
  
  gtf_filt <- gtf_raw[gtf_rows,]
  return(gtf_filt)
}

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

fit_ld <- function(marey,ld_file){
  marey_narrow <- marey %>% dplyr::select(SNPname,mb,homolog)
  
  ld.bychr.raw <- fread(ld_file)
  
  ld.bychr <- ld.bychr.raw %>% dplyr::select(SNP_A,SNP_B,R2) %>% left_join(marey_narrow,by=c("SNP_A"="SNPname")) %>% left_join(marey_narrow,by=c("SNP_B"="SNPname"))  %>% filter(homolog.x == homolog.y)
  
  ld.bychr$mb_distance <- abs(ld.bychr$mb.x - ld.bychr$mb.y)
  ld.bychr$R2 <- as.numeric(ld.bychr$R2)
  
  ld.bychr$chrType[ld.bychr$homolog.x != "Z"] <- "A"
  ld.bychr$chrType[ld.bychr$homolog.x == "Z"] <- "Z"
  
  fit_hyperbola <- function(chrom,ld.bychr){
    fit_tmp <- lm(na.omit(ld.bychr$R2[ ld.bychr$chrType == chrom]) ~ I(1/na.omit(ld.bychr$mb_distance[ ld.bychr$chrType == chrom])))
    
    ld_hyperbola <- function(x){
      y = fit_tmp$coefficients[2]/x + fit_tmp$coefficients[1]
    }
    return(ld_hyperbola)
  }
  
  ld_hyperbola_A <- fit_hyperbola(chrom = "A",ld.bychr)
  ld_hyperbola_Z <- fit_hyperbola(chrom = "Z",ld.bychr)
  
  pdf(paste("fig_S2.pdf",sep=''),width=5.5,height=4)
  print(
    ggplot(ld.bychr ,aes(x=mb_distance,y=R2,color=as.factor(chrType))) +  
      geom_point(aes(alpha=ifelse(chrType == 'Z',1,0.01),shape=ifelse(chrType == 'Z',"Z","A")),size=0.5) + 
      theme_bw()  + ylim(0,1)+ xlim(0,0.1)+ 
      geom_function(fun = ld_hyperbola_A, colour = "red") +
      geom_function(fun = ld_hyperbola_Z, colour = "blue") + 
      scale_color_manual(values=c("red","blue")) + 
      guides(alpha='none',shape='none')+
      scale_shape_manual(values=c(20,19)) +
      scale_alpha(range = c(0.01,0.5)) +
      labs(x="Pairwise Distance (mb)",
           y=expression(paste("Linkage Disequilibrium (",R^2,")")),color="Chrom. Type") 
  )
  dev.off()
  
}
