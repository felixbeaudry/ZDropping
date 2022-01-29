#modified from Windowing_genome_content.R by Joanna.Rifkin on Nov 25 2020

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

#clustering packages
library(clustertend)
library(cluster)
library(fpc)
library(corrplot)
library(pvclust)
library(factoextra)
library(fpc)
library(NbClust)
library(Rtsne)

source("~/Google Drive/Research/Scripts/linked-selection-master/linkedsel_functions.R")
est_rec_piece_func  <- function(x) {est_rec_piece(x)}

####functions####
'%ni%' <- function(x,y)!('%in%'(x,y))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

assemblySizes <- function(sizesFile){  
  fread(sizesFile)
  
  names(sizes) <- c("scaff","bp")
  
  cat("genome size is: ",sum(sizes$bp),'bp\n')
  
  print(ggplot(sizes,aes(x=bp/1000000)) + geom_histogram(bins=30))
  return(sizes)
}

windowedMap <- function(window_size=1000000,sizes){
  
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
  
  for(chrom in unique(marey_zerod$crowmarey_zerod$homolog != "Unplaced")){
    marey_temp <- marey_zerod[marey_zerod$homolog == chrom,]
    if(length(unique(marey_temp$cm))>4){
      if(cor.test(marey_temp$cm,marey_temp$mb)[4] < 0){
        marey_zerod$cm_flip[marey_zerod$crow == chrom] <- 
          max(marey_zerod$cm[marey_zerod$crow == chrom]) - marey_zerod$cm[marey_zerod$crow == chrom]
      }
    }
  }
  
  marey_zerod <- marey_zerod %>% filter(!is.na(cm_flip) & !is.na(mb) & homolog %ni% c("W","Unplaced"))
  
  marey_max <-  marey_zerod %>% group_by(homolog) %>% summarise(max_cm = max(cm_flip), max_mb=max(mb))
  marey_max$cmmb <- marey_max$max_cm/marey_max$max_mb
  marey_max<-  marey_max[complete.cases(marey_max),]
  marey_max$max_log_mb <- log(marey_max$max_mb)
  
  
  fullMap_plot <- 
    ggplot(marey_zerod ,aes(x=mb,y=cm_flip)) +
    
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
ldhatFit <- function(ldhatFile){
  load(ldhatFile)
  
  #subsample
  LD.list.ordered_sub <- LD.list.ordered[[1]][seq(1, nrow(LD.list.ordered[[1]]), 100), ]
  
  #convert into df
  for(i in c(2:length(LD.list.ordered))){
    LD.list.ordered_sub <- rbind.data.frame(LD.list.ordered_sub,
                                            LD.list.ordered[[i]][seq(1, nrow(LD.list.ordered[[i]]), 100), ]
    )
  }
  #LD.list.ordered_sub <- separate(LD.list.ordered_sub, CHROM, c("dovetail","crow"),   sep = "ch", remove = TRUE, convert = FALSE, extra = "merge", fill = "right")
  
  LD.list.ordered_sub$mb <- LD.list.ordered_sub$POS / 1000000
  LD.list.ordered_sub$Cum_rho_adj <- LD.list.ordered_sub$Cum_rho/5000 #approx Ne?
  return(LD.list.ordered_sub)
}
compareRHO2MB <- function(windowedreC_join){
  windowedreC_join <- windowedreC_join[complete.cases(windowedreC_join),]
  
  #filter
  windowedreC_join <- windowedreC_join %>% filter(rho > 2 | link < 6 | xor(rho < 2 , link > 6))
  
  #calculate correlations for each window
  rec_corr <- foreach(chrom=unique(windowedreC_join$LG),.combine=rbind) %do% {
    windowedreC_tmp <- windowedreC_join %>% filter(LG == chrom)
    rec_corr_tmp <- c(chrom,cor(windowedreC_tmp$rho,windowedreC_tmp$link,method = "spearman"))
    rec_corr_tmp
  }
  rec_corr <- as.data.frame(rec_corr)
  names(rec_corr) <- c("chrom","corr")
  rec_corr$corr <- as.numeric(as.character(rec_corr$corr))
  
  windowedreC_join <- left_join(windowedreC_join,rec_corr,by=c("LG"="chrom")) %>% left_join(sizes,by=c("LG"="scaf"))
  windowedreC_join$chr <- factor(windowedreC_join$LG, levels=c("1"  ,"1A","2" , "3"  ,"4" ,"4A", "5","Z",  "6",  "7",  "8" , "9" ,"10" ,"11", "12" ,"13", "14", "15" ,"16", "17" ,"18", "19", "20" ,"21","22", "23", "24", "25", "26", "27", "28"))
  
  print(
    ggplot(windowedreC_join,aes(x=link,y=rho)) + 
      geom_point(alpha=0.75,size=0.75) + 
      stat_smooth(method="lm",formula = 'y~x',se=FALSE,aes(color=corr))+  
      facet_wrap(. ~chr,scales="free") + theme_bw() +
      labs(x="Linkage (cM/mb)",y="LD (rho)",color="Corr. (Sp.Rank)") + 
      scale_color_gradient(
        low = "white",
        high = "black",
        na.value = "gray"
      ) 
  )
}

makeSexMarey <- function(sexlinkagefile,marey_zerod){
  sex_map_raw <-  fread(sexlinkagefile)[,-10]
  names(sex_map_raw) <- c("lg","markerNum","marker","male","female","male_r","male_r_lnl","female_r","female_r_lnl")
  
  #remove unlikely markers
  sex_map <- sex_map_raw %>% filter(marker %in% genmap_filt$marker) 
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
  
  sex_map_melt <- melt(sex_map[,c(1,3,10:11)],id.vars = c("marker","lg"),verbose=FALSE)
  
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
    #cat(chr," ",max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable == "male"],na.rm=T)," ")
    totalLength_m = max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable == "male"],na.rm=T) + totalLength_m
    totalLength_f = max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable == "female"],na.rm=T) + totalLength_f
    if(max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable == "male"],na.rm=T) > max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable == "female"],na.rm=T)){cat(chr,":male ")}
    if(max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable == "male"],na.rm=T) < max(sex_map_melt$value[sex_map_melt$lg == chr & sex_map_melt$variable == "female"],na.rm=T)){cat(chr,":female ")}
    
  }
  
  ##with mb position
  sex_marey <- left_join(sex_map_melt,marey_zerod[,c(1,13,12)],by=c("marker"="SNPname"))
  
  ##flip chromosome SNP mb positionss so they follow the marey map
  sex_marey_zerod <-  sex_marey[order(sex_marey$lg,sex_marey$variable,sex_marey$value),]
  sex_marey_zerod$cm_flip <- sex_marey_zerod$value
  
  for(chrom in unique(sex_marey_zerod$lg)){
    # cat(chrom)
    marey_temp <- sex_marey_zerod[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name == "male",] #chose male so flip can work on Z too
    
    if(cor.test(marey_temp$value,marey_temp$mb)[4] < 0){
      sex_marey_zerod$cm_flip[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "male"] <- 
        max(sex_marey_zerod$value[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "male"]) - sex_marey_zerod$value[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable  == "male"]
      
      sex_marey_zerod$cm_flip[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "female"] <- 
        max(sex_marey_zerod$value[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable_name  == "female"]) - sex_marey_zerod$value[sex_marey_zerod$lg == chrom & sex_marey_zerod$variable  == "female"]
    }
  }
  
  sex_marey_zerod <- sex_marey_zerod[complete.cases(sex_marey_zerod),]
  
  cat("Plot 1/2: Male to Female Recombination Rate")
  print(sex_dist_plot)
  cat("Plot 2/2: Male and Female Genetic Distances")
  print(sex_map_plot)
  return(sex_marey_zerod)
}

makeLinkageSDI <- function(rec_windows_male,rec_windows_female){
  names(rec_windows_male)[3] <- "male"
  names(rec_windows_female)[3] <- "female"
  
  rec_windows_sex <- left_join(rec_windows_male,rec_windows_female)
  
  rec_windows_sex$male <- as.numeric(rec_windows_sex$male)
  rec_windows_sex$female <- as.numeric(rec_windows_sex$female)
  
  rec_windows_sex$SDI <- ifelse(rec_windows_sex$male >rec_windows_sex$female,-1*((rec_windows_sex$male/rec_windows_sex$female) - 1), 1*((rec_windows_sex$female/rec_windows_sex$male) - 1))
  return(rec_windows_sex)
}

calcGapsNBases <- function(bybaseFile,homology_conversion,  gap_window= 1000000){
  
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
  
  bybase_melt <- reshape2::melt(bybase_filt,id.vars = c("homolog","mb"),verbose=FALSE)
  bybase_melt$type <- factor(bybase_melt$variable, levels=c( "GC", "AT", "N_frac"))
  return(bybase_melt)
}

makeCpG <- function(cpgFile,homology_conversion){
  
  
  cpg_tally <- fread(cpgFile) %>% group_by(chrom) %>% tally() %>% left_join(homology_conversion,by=c("chrom"="scaff"))
  
  cpg_tally$cpgMb <- cpg_tally$n / (cpg_tally$bp/1000000)
  
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
makeRepeatWindows <- function(repeatFile, window_size_TE = 1000000,df_map,homology_conversion){
  
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
  
  
  All_repeats_win <-  
    All_repeats_mut %>% group_by(scaff,family)  %>% 
    mutate(position_window=LG_lower_start_adj%/%as.numeric(window_size_TE)) %>% 
    group_by(scaff,position_window,family) %>%  
    add_tally() %>% select_if(., is.numeric) %>% 
    summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) %>% 
    mutate(proportion_of_window=featuresize_sum/as.numeric(window_size_TE)) %>% 
    mutate(window_start=(position_window*(as.numeric(window_size_TE)+1)) )
  
  All_repeats_win$mb <- round(as.numeric(All_repeats_win$window_start) / 1000000)
  
  ##add in ZERO windows for each family!
  #make data.frame with full length of windows for each chrom/scaff -- window map from Hi-C ??
  
  
  
  
  df_map_TE <- 
    foreach(TEfam=unique(All_repeats_win$family),.combine=rbind.data.frame) %do% {
      cbind.data.frame(df_map,"TEfam"=TEfam )
    }
  
  All_repeats_zeros <- left_join(df_map_TE,All_repeats_win,by=c("V1"="scaff","mb"="mb","TEfam"="family"))
  
  All_repeats_zeros$proportion_of_window[is.na(All_repeats_zeros$proportion_of_window)] <- 0
  
  All_repeats_zeros <- left_join(All_repeats_zeros,homology_conversion,by=c("V1"="scaff"))
  
  
  
  print(
    # ggplot(All_repeats_zeros %>% filter(!is.na(homolog ) & homolog %ni% c("16","W","29") & proportion_of_window <= 1),aes(x=mb,y=proportion_of_window,fill=TEfam)) + 
    ggplot(All_repeats_win ,aes(x=mb,y=proportion_of_window,fill=family)) + 
      
      geom_area(position="stack") + #ylim(0,0.2)+
      # scale_x_continuous(breaks=c(50,100))+
      
      facet_wrap(. ~ scaff, scales = 'free',strip.position = "bottom") +
      labs(x="", y=" ",fill=" ") +
      theme_bw(base_size=8) + 
      #guides(fill="none")+
      theme(strip.background =element_rect(fill="white"),axis.title.x=element_blank()) + #ylim(0,1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      theme(panel.spacing = unit(-0.1, "cm"),
            panel.border = element_rect(color = "white", fill = NA, size = 0))  +
      #  coord_cartesian(ylim = c(0, 0.3))+
      scale_fill_manual(values=c("#ffd700",
                                 "#ffb14e",
                                 "#fa8775",
                                 "#ea5f94",
                                 "#cd34b5",
                                 "#9d02d7",
                                 "#0000ff","grey"))
  )
  return(All_repeats_zeros)
}
####sizes####
setwd('~/Google Drive/Research/Data2/fsj/genome')

window_size = 1000000

sizes <- assemblySizes(sizesFile='FSJ.chrom.sizes') #summarize sizes

df_map <- windowedMap(sizes=sizes) #split each scaffold into wins (of window_size)

size_homology <- homologyMap(conversionFile='FSJV3_convert.txt',size=sizes) #Chicken-based homology

####Recombination####

genmap <- makeGeneticMap(genmapFile='crimap/all.fixed.txt')
frameworkSNPs <- fread('frameworkSNPs.list')
frameworkSNPs$framework <- 1

genmap_frame <- left_join(frameworkSNPs,genmap,by=c("V1"="lg","V2"="raworder"))
names(genmap_frame)[c(1,2)] <- names(genmap)[c(1,2)]

genmap_boot <- genmap_frame
genmap_boot <- genmap


w_size_cm = 10
#nloci<- 5000 #set sim loci number relative to window size; regular sim nloci / genome size in mb

#loop across scaffolds making windows of SNPs
win_global = 0

genmap_boot$bootstrap <- NA

for(LG in unique(genmap_boot$lg)){
  
  size_tmp <- na.omit(genmap_boot$cm[genmap_boot$lg == LG])
  
  #loop across windows
  for(win in seq(from=0,to=max(genmap_boot$cm[genmap_boot$lg == LG],na.rm=T),by = w_size_cm)){
    if(win+w_size_cm<max(size_tmp)){
      cat(LG," ", ((win/w_size_cm) + win_global)," ",win," ",size_tmp,"\n")
      genmap_boot$bootstrap[genmap_boot$cm > win & genmap_boot$cm <= (win+w_size_cm) & genmap_boot$lg == LG] <- (win/w_size_cm) + win_global
    } #else{ #skip windows at the end of scaffolds with too few SNPs
     # cat(LG," ", ((win/w_size_cm) + win_global)," ",win," ",size_tmp," skipped\n")
      
   # }
  }
  win_global <- (win/w_size_cm) + win_global #need to increase window number with each loop across chromosomes
  
}



genmap_boot_tally <- genmap_boot %>% group_by(bootstrap,lg) %>% tally()
genmap_boot_tally %>% group_by(lg) %>% summarise(ave_SNP = mean(n))


ggplot(genmap_boot_tally %>% filter(!is.na(bootstrap)),aes(x=n,fill=as.factor(lg))) + geom_histogram()

##

marey <- makeMareyMap(chipFile='FSJV3.beadChip.r060122.sam.pos',genmap=genmap,size_homology=size_homology)

rec_windows_cm <- get_rec_windows(CMMB = marey[,c(15,1,18,13)]) #columns = "chr"    "marker" "cm"     "mb"    

ldhat_processed <- ldhatFit(ldhatFile='percontig_LDhat_rho.Rdata')

rec_windows_rho <- get_rec_windows(CMMB = ldhat_processed[,c(1,3,10,9)]) # cols: "chr", "marker", "cm", "mb"

#compare linkage and rho
names(rec_windows_rho)[3] <- "rho"
names(rec_windows_cm)[3] <- "link"

windowedreC_join <- left_join(rec_windows_cm,rec_windows_rho)

compareRHO2MB(windowedreC_join=windowedreC_join)

#sex & linkage
sex_marey_zerod <- makeSexMarey(sexlinkagefile='crimap/all.sex.fixed.txt',marey_zerod=marey)

rec_windows_male <- get_rec_windows(CMMB = sex_marey_zerod[sex_marey_zerod$variable == "male_adj" ,c(2,1,10,9),]) # cols: "lg", "marker", "cm", "mb"
rec_windows_female <- get_rec_windows(CMMB = sex_marey_zerod[sex_marey_zerod$variable == "female_adj" ,c(2,1,10,9),]) # cols: "lg", "marker", "cm", "mb"

rec_windows_sex <- makeLinkageSDI(rec_windows_male=rec_windows_male,rec_windows_female=rec_windows_female)

#sequence content
bybase <- calcGapsNBases(bybaseFile='gapsIN.FSJ.v3.win1mb.txt',homology_conversion=size_homology)

cpg <- makeCpG(cpgFile='all.FSJ3.cpg',homology_conversion=size_homology) ##cpg islands##

repeats_windowed <- makeRepeatWindows(repeatFile='FSJ.hifiasm.scaff.fasta.out',df_map=df_map,homology_conversion=size_homology)


####HiC####
hic_raw <- fread('FSJV3.hic.1mb.raw.txt')

hic_inter <- hic_raw %>% group_by(V1,V2) %>% 
  summarise(mean_connect = mean(V5)) %>% 
  left_join(size_homology,by=c("V1"="scaff")) %>% left_join(size_homology,by=c("V2"="scaff")) 

hic_filt <- hic_inter %>% 
  filter(homolog.x != homolog.y & homolog.x %ni% c("16","W") & homolog.y %ni% c("16","W") & !is.na(homolog.y) & !is.na(homolog.x))

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


#see reformat_hic.R
load('hic_m.1mb.raw.rdata')

hic_pca <- PCA(hic_m_n2,graph = F) #should this be done within chromosomes only?
hic_pca_coord <- cbind.data.frame(hic_df_map,as.data.frame(hic_pca$var$coord))
hic_pca_coord$pos <- as.numeric(hic_pca_coord$V2) / 1000000

ggplot(hic_pca_coord,aes(x=pos,y=Dim.1)) + 
  geom_point(alpha=0.2) + geom_line() +  facet_wrap(Chr ~ .,scales="free") + 
  theme_bw() 

##as eigenvectors
#hic_eigen <- eigen(hic_m_n2)
#hic_eigen$vectors

#reformat
hic_m_melt <- melt(hic_m_n2)

#add/join position information
hic_m_map <- left_join(hic_m_melt,hic_df_map,by=c("Var1"="mat_pos")) %>% left_join(hic_df_map,by=c("Var2"="mat_pos"))

hic_m_map <- hic_m_map[c(7,5,12,10,3)]
names(hic_m_map) <- c("chr_1","pos_1","chr_2","pos_2","links_norm")

hic_m_map$connectType[hic_m_map$chr_1 == hic_m_map$chr_2 ] <- "intra"
hic_m_map$connectType[hic_m_map$chr_1 != hic_m_map$chr_2 ] <- "inter"

hic_byChrom <- 
  hic_m_map %>% filter(chr_1 %in% sizes$chrom[sizes$chromType %in% c("micro","macro","int")] & chr_2 %in% sizes$chrom[sizes$chromType %in% c("micro","macro","int")]) %>% group_by(chr_1,chr_2,connectType) %>% summarise(mean_connect = mean(links_norm))

hic_byChrom$Chr_1 <- factor(hic_byChrom$chr_1, levels=size_order$chrom)
hic_byChrom$Chr_2 <- factor(hic_byChrom$chr_2, levels=size_order$chrom[order(size_order$size)])

  ggplot(hic_byChrom %>% filter(connectType == "inter"), aes(x=Chr_1, y=Chr_2, fill= log(mean_connect))) + 
  geom_tile(color="lightgray") + theme_bw(base_size=15) + 
  scale_fill_gradient2(
    low = "cornflowerblue",
    mid = "white",
    high = "indianred",
    midpoint = median(log(hic_byChrom$mean_connect[hic_byChrom$connectType == "inter"])),
    na.value = "gray"
  ) + labs(y="",fill="log norm. Hi-C links") + #guides(fill=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank())   +
  theme(legend.position="bottom") + theme(aspect.ratio=1) +
    scale_x_discrete(guide = guide_axis(n.dodge=2))

#along chromosomes
hic_bypos_connect <- hic_m_map %>% group_by(chr_1,pos_1,connectType) %>% summarise(mean_connect = mean(links_norm))

hic_bypos_connect$mb <- as.numeric(hic_bypos_connect$pos_1) / 1000000

hic_bypos_connect <- left_join(hic_bypos_connect,sizes,by=c("chr_1"="chrom"))

hic_bypos_connect$chrom <- factor(hic_bypos_connect$chr_1, levels=c("1"  ,"1A","2" , "3"  ,"4" ,"4A", "5","Z",  "6",  "7",  "8" , "9" ,"10" ,"11", "12" ,"13", "14", "15" ,"16", "17" ,"18", "19", "20" ,"21","22", "23", "24", "25", "26", "27", "28"))

hic_bypos_connect_i <- hic_bypos_connect[hic_bypos_connect$connectType == "inter",]


  ggplot(hic_bypos_connect_i[hic_bypos_connect_i$chromType %in% c("macro"),],aes(x=mb,y=mean_connect,fill=connectType)) + 

  geom_area(position="stack") + #ylim(0,0.2)+
  scale_x_continuous(breaks=c(50,100))+
  
  facet_grid(. ~ chrom, scales = 'free_x', space = 'free_x', switch="x") +
  labs(x="", y=" ",fill=" ") +
  theme_bw(base_size=15) + 
  guides(fill=FALSE)+
  theme(strip.background =element_rect(fill="white"),axis.title.x=element_blank()) + #xlim(120,125) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(panel.spacing = unit(-0.1, "cm"),
        panel.border = element_rect(color = "white", fill = NA, size = 0)) 


####by chrom analyses####
sizes_chrom <- sizes %>% filter(chrom %ni% c("29","U","W","16","31_LG33_ZFLG2","LG32_ZFLGE22","LG34_ZFUn","LG36","LG37","LG43","30","31","32","34","35","37"))
sizes_chrom <- sizes_chrom[,c(3,2)]

repeats_grouped <- 
  All_repeats_win %>% filter(family %in% c("DNA_transposon","house_RNAs","LINE","Low_complexity","LTR","Non-LTR","Satellite","Simple_repeat"),chrom %in% sizes$chrom,chrom %ni% c("U","W","16","31_LG33_ZFLG2","LG32_ZFLGE22","LG34_ZFUn","LG36","LG37","LG43","29","30","31","32","34","35","37"),) %>% group_by(chrom,family) %>% summarize(win_mean=mean(proportion_of_window)) %>% dcast(chrom~family)

#coding_count_mb_chrom <- coding_count_mb %>%   group_by(chrom) %>% tally()
#genes_byChrom <- left_join(
#  stringtie_byChrom[,c(7,8)],
#  coding_count_mb_chrom,by=c("chrom"="chrom"))
names(cpg_tally)[3] <- "CpG"
#names(genes_byChrom)[3] <- "CDS"

hic_inter_byChrom <- 
  hic_byChrom %>% filter(connectType == "inter") %>% group_by(chr_1) %>% summarise(hic=mean(mean_connect))

table_tmp1 <- left_join(sizes_chrom,
                        unique(bybase[,c(1,3,4)] %>% group_by(chrom) %>% summarise(ave_GC = ave(GC))),
                        by=c("chrom"="chrom")) %>%
  left_join(cpg_tally[,c(1,4)], by=c("chrom"="chrom")) %>%
  left_join(repeats_grouped,by=c("chrom"="chrom")) %>%
#  left_join(genes_byChrom,by=c("scaf"="chrom")) %>%
  left_join(
    marey_zerod %>% group_by(crow) %>% summarise(map_length=max(cm)),by=c("chrom"="crow"))

allOfit <-  left_join(table_tmp1,sex_max[-c(22,24),c(4,5)],by=c("chrom"="crow")) %>%
  left_join(hic_inter_byChrom,by=c("chrom"="chr_1"))

allOfit_norm <- allOfit
#allOfit_norm[,c(5:13,15:16)] <- allOfit_norm[,c(5:13,15:16)] / allOfit_norm$chrom_length
allOfit_norm$map_length <- allOfit_norm$map_length / allOfit_norm$chrom_length

allOfit_scale <- scale(allOfit_norm[,-1])
row.names(allOfit_scale) <- allOfit_norm$chrom
colnames(allOfit_scale) <- 
  c("Length",
    #"Gaps",
    "GC",
    "CpG",
    "DNA transposons",
    "house RNAs",
    "LINEs",
  #  "Low complex.",
  "Non-LTR",
    "LTR",
    "Satellites",
    "Simple repeats",
  #  "SINE",
    #"Transcription",
   # "Coding Seq",
    "Recombination",
    "Sex*Rec.",
    "Hi-C links")
    

save(allOfit_norm,file="allOfit_norm_V3.rdata")

#load("allOfit_norm.rdata")
##
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
# matrix of the p-value of the correlation
p.mat <- cor.mtest(allOfit_scale)
allOfit_scale_M<-cor(allOfit_scale)

par(mfrow=c(1,1))
corrplot(allOfit_scale_M, method="color",  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE ,tl.cex = 0.9 ,number.cex = 0.6
)

##tsne
#Labels<-allOfit_norm$scaf
#allOfit_norm$scaf<-as.factor(allOfit_norm$scaf)
## for plotting
#colors = rainbow(length(unique(allOfit_norm$scaf)))
#names(colors) = unique(allOfit_norm$label)

## Executing the algorithm on curated data
#tsne <- Rtsne(allOfit_scale, dims = 2, perplexity=3, verbose=TRUE, max_iter = 5000)

## Plotting
#plot(tsne$Y, t='n', main="tsne")
#points(tsne$Y, col=colors[allOfit_norm$chrom_length])
#text(tsne$Y, labels=allOfit_norm$scaf, col=colors[allOfit_norm$scaf])
#text(tsne$Y, labels=allOfit_norm$chrom)

#tsne_tmp <- cbind.data.frame("dim.1"=tsne$Y[,1],"dim.2"=tsne$Y[,2],"chrom"=allOfit_norm$chrom,"length"=allOfit$chrom_length)

#tSNE_plot <- 
  
#ggplot() +
#  geom_point(data=tsne_tmp,aes(x=dim.1,y=dim.2,color=length/1000000),size=8) + 
#  theme_bw() +
#  geom_text(data=tsne_tmp,aes(x=dim.1,y=dim.2,label=chrom),color="white") +
#  guides(color=FALSE)+
#  scale_colour_gradient(low = "skyblue", high = "darkgreen", na.value = NA) +
#  labs(color="Size (mb)",size="log(Size)",x="",y="",title="t-SNE") 

##PCA
allOfit_pca <- PCA(allOfit_scale,graph = F)

allOfit_pca_labs <- cbind.data.frame("chrom"=allOfit_norm$chrom,allOfit_pca$ind$coord,"chrom_length"=allOfit_norm$chrom_length)
allOfit_pca_contribs<- cbind.data.frame("Var"=row.names(allOfit_pca$var$coord),allOfit_pca$var$coord)

allOfit_pca$var$coord
allOfit_pca$var$contrib

allOfit_pca_contribs$Varf <- factor(allOfit_pca_contribs$Var,levels=allOfit_pca_contribs$Var[order(-allOfit_pca_contribs$Dim.1)])

plot_grid(
  
  ggplot(allOfit_pca_contribs,aes(x=Varf,y=Dim.1)) + geom_bar(stat="identity") + theme_bw(base_size = 14) + labs(x="",y="PC1 % Variance") +
    scale_x_discrete( guide = guide_axis(n.dodge=2)) ,
  ggplot(allOfit_pca_contribs,aes(x=Varf,y=Dim.2)) + geom_bar(stat="identity") + theme_bw(base_size = 14) +  labs(x="",y="PC2 % Variance") +
    scale_x_discrete( guide = guide_axis(n.dodge=2)) ,
  labels = c('', ''), ncol = 1, align = 'hv',axis='tbrl')

#PCA_plot <- 
ggplot() +
  geom_point(data=allOfit_pca_labs,aes(x=Dim.1,y=Dim.2,color=chrom_length/1000000),size=8) + 
  theme_bw() +
 geom_text(data=allOfit_pca_labs,aes(x=Dim.1,y=Dim.2,label=chrom),color="white") +
 #  guides(color=FALSE)+
  scale_colour_gradient(low = "skyblue", high = "darkgreen", na.value = NA) +
  labs(color="Size (mb)",size="log(Size)",title="PCA",
       x=paste("PC1 (",round(allOfit_pca$eig[1,2],digits=2),"%)",sep=""),y=paste("PC2 (",round(allOfit_pca$eig[2,2],digits=2),"%)",sep="")) +
  coord_sf()
  
  ggrepel::geom_text_repel(
    min.segment.length = 0.5,
    data = allOfit_pca_contribs,
    mapping = aes(x=Dim.1*5,y=Dim.2*5,label = Var),
    position = position_dodge(1.5)
  )
  
#  plot_grid(tSNE_plot,  PCA_plot,  labels = c('', ''), ncol = 2, align = 'hv',axis='tbrl')
  
  ggplot() +
  #  geom_point(data=allOfit_pca_labs,aes(x=Dim.1,y=Dim.2,color=chrom_length/1000000),size=1) + 
    theme_bw() +
#    geom_text(data=allOfit_pca_labs,aes(x=Dim.1,y=Dim.2,label=scaf),color="white") +
    guides(color=FALSE)+
    scale_colour_gradient(low = "skyblue", high = "darkgreen", na.value = NA) +
    labs(color="Size (mb)",size="log(Size)",title="PCA",
         x=paste("PC1 (",round(allOfit_pca$eig[1,2],digits=2),"%)",sep=""),y=paste("PC2 (",round(allOfit_pca$eig[2,2],digits=2),"%)",sep="")) +
  
  ggrepel::geom_text_repel(
    min.segment.length = 0.5,
    data = allOfit_pca_contribs,
    mapping = aes(x=Dim.1,y=Dim.2,label = Var),
    position = position_dodge(1.5)
  ) + coord_sf()
  
##k-means
#optimal k
wss <- (nrow(allOfit_scale)-1)*sum(apply(allOfit_scale,2,var))
for (i in 2:27) wss[i] <- sum(kmeans(allOfit_scale,
                                     centers=i)$withinss)
plot(1:27, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#pam
v <- rep(0, 10)
  for (i in 2:10) {
    clust <- pam(allOfit_scale, i)
    ss <- cluster::silhouette(clust$cluster, stats::dist(allOfit_scale))
    v[i] <- mean(ss[, 3])
    # .get_ave_sil_width(, )
  }

pam_sil <- cbind.data.frame(v,"k"=seq(1,length(v),1))
ggplot(pam_sil,aes(x=k,y=v))+geom_point() + geom_line() +
  geom_vline(xintercept = 2,linetype="longdash",color="gray",alpha=0.75)+
  theme_bw(base_size=15) +labs(y="Average silhouette width",x="Number of cluster k")

pam(allOfit_scale,method="silhouette")

nb <- NbClust(allOfit_scale, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans")

#

km.res2 <- kmeans(allOfit_scale, centers=2,iter.max = 20)
km.res3 <- kmeans(allOfit_scale, centers=3,iter.max = 20)


#clusplot(allOfit_scale, km.res2$cluster, color=TRUE, shade=FALSE, labels=2, lines=0,)

#fviz_cluster(list(data = allOfit_scale, cluster = km.res2$cluster),
#             ellipse.type = "norm", geom = "text", stand = FALSE, 
#             palette = "jco", ggtheme = theme_classic(),repel=TRUE)

allOfit_pca_labs$cluster2 <- km.res2$cluster
allOfit_pca_labs$cluster3 <- km.res3$cluster

##
ggplot(data=allOfit_pca_labs,aes(x=Dim.1,y=Dim.2,color=as.factor(cluster3))) +
  geom_point(size=8) + 
  theme_bw() +
  geom_text(aes(label=chrom),color="white") +
  guides(color=FALSE)+
 
#  stat_ellipse(type = "t", linetype = 3) +
 # stat_ellipse(type = "norm", linetype = 2) + 
  coord_fixed() +
  stat_ellipse(type = "euclid", level = 3) +  #circle
  
  scale_colour_manual(values = c( "skyblue",  "darkgreen" ,"#cd34b5")) +
  labs(color="Size (mb)",size="log(Size)",title="PCA",
       x=paste("PC1 (",round(allOfit_pca$eig[1,2],digits=2),"%)",sep=""),y=paste("PC2 (",round(allOfit_pca$eig[2,2],digits=2),"%)",sep=""))# +

ggrepel::geom_text_repel(
  min.segment.length = 0.5,
  data = allOfit_pca_contribs,
  mapping = aes(x=Dim.1*5,y=Dim.2*5,label = Var),
  position = position_dodge(1.5)
)



#pamk(allOfit_scale)
#pam.res <- pam(allOfit_scale, 2)

#fviz_cluster(pam.res,
#             palette = c("#00AFBB", "#FC4E07"), # color palette
#             ellipse.type = "t", # Concentration ellipse repel = TRUE, # Avoid label overplotting (slow) ggtheme = theme_classic()
#)+ theme_classic()

# Centroid Plot against 1st 2 discriminant functions
#plotcluster(allOfit_scale, km.res2$cluster)


## hierarchical clustering

d <- dist(allOfit_scale, method = "euclidean") # distance matrix
fit_hclust <- hclust(d, method="ward") #hclust is Agglomerative clustering, “bottom-up” manner, good at small groups

plot(fit_hclust) # display dendogram
groups <- cutree(fit_hclust, k=2)
rect.hclust(fit_hclust, k=2, border="red")

res.diana <- diana(x = allOfit_scale, # data matrix
                   stand = TRUE, # standardize the data
                   metric = "euclidean"# metric for distance matrix 
)

#fviz_dend(res.diana, cex = 0.6, k = 2)       



fit_pvclust <- pvclust(t(allOfit_scale), method.hclust="ward",
               method.dist="euclidean")
plot(fit_pvclust,c("au"),print.num=F,col.pv=c(au=4)) # dendogram with p values
pvrect(fit_pvclust, alpha=.95) # add rectangles around groups highly supported by the data


library(mclust)



v <- rep(1, 8)
for (i in 1:8) {
  fit_mclust <- Mclust(allOfit_scale,G=i)
  v[i] <- max(fit_mclust$BIC,na.rm = T)
  # .get_ave_sil_width(, )
}
mclust_BIC <- cbind.data.frame(v,"k"=seq(1,length(v),1))
ggplot(mclust_BIC,aes(x=k,y=v))+geom_point() + geom_line() +
  geom_vline(xintercept = 3,linetype="longdash",color="gray",alpha=0.75)+
  theme_bw(base_size=15) +labs(y="Bayesian Information Criterion (BIC)",x="Number of groups G")

fit_mclust_G3 <- Mclust(allOfit_scale,G=3)
allOfit_pca_labs$mcluster3 <- fit_mclust_G3$classification  
#fit_mclust_G4 <- Mclust(allOfit_scale,G=4)
#allOfit_pca_labs$mcluster4 <- fit_mclust_G4$classification  

##
ggplot(data=allOfit_pca_labs,aes(x=Dim.1,y=Dim.2,color=as.factor(mcluster3))) +
  geom_point(size=8) + 
  theme_bw() +
  geom_text(aes(label=chrom),color="white") +
  guides(color=FALSE)+
  
  #  stat_ellipse(type = "t", linetype = 3) +
  stat_ellipse(type = "norm") + coord_fixed() +
#  stat_ellipse(type = "euclid", linetype = 2, level = 4) +  #circle
  
  #  scale_colour_manual(values = c( "skyblue",  "darkgreen")) +
  labs(color="Size (mb)",size="log(Size)",title="PCA G=3",
       x=paste("PC1 (",round(allOfit_pca$eig[1,2],digits=2),"%)",sep=""),y=paste("PC2 (",round(allOfit_pca$eig[2,2],digits=2),"%)",sep=""))# +

ggrepel::geom_text_repel(
  min.segment.length = 0.5,
  data = allOfit_pca_contribs,
  mapping = aes(x=Dim.1*5,y=Dim.2*5,label = Var),
  position = position_dodge(1.5)
)





####by window####

#
All_repeats_win_dcast <- 
  All_repeats_scafs[,c(1,3,30,28)] %>%
dcast(chr*mb~family,value.var= "proportion_of_window")

All_repeats_win_dcast[is.na(All_repeats_win_dcast)] <- 0
All_repeats_win_dcast$mb <- round(All_repeats_win_dcast$mb,digits = 0)



bybase_win <- bybase[,c(1:5)]
cpg_win <- cpg[,c(2,4)]
cpg_win$mb <- round(cpg_win$From/1000000,0)
cpg_wint <- cpg_win %>% group_by(scaf,mb) %>% tally()
names(cpg_wint)[3] <- "cpgs"


coding_count_win <- coding_count_mb_gene[,c(1:3)]
names(coding_count_win)[3] <- "genes"

windowedreC_ave_win <- windowedreC_join[,c(1,4,6,8)]

windowedreC_sex_win <- windowedreC_sex[,c(11,2,3,4,5)]

hic_win <- hic_bypos_connect[,c(10,5,3,4)] %>%
  dcast(chrom*mb~connectType,value.var= "mean_connect")

all_the_wins <- left_join(hic_win,windowedreC_sex_win,by=c("chrom"="Chrom","mb"="mb_s")) %>%
  left_join(windowedreC_ave_win,by=c("chrom"="chr","mb"="mb_s"))  %>%
left_join(coding_count_win,by=c("chrom"="chrom","mb"="mb_r"))%>%
  left_join(cpg_wint,by=c("chrom"="scaf","mb"="mb"))%>%
  left_join(bybase_win,by=c("chrom"="chrom","mb"="mb"))%>%
  left_join(All_repeats_win_dcast,by=c("chrom"="chr","mb"="mb"))

all_the_wins$label <- paste("chr",all_the_wins$chrom,"mb",all_the_wins$mb,sep="")
  
all_the_wins_comp <- all_the_wins[complete.cases(all_the_wins),]
all_the_wins_comp <- all_the_wins_comp[,-c(14,20,24)]

all_the_wins_scale <- scale(all_the_wins_comp[,-c(1,2,22)])
#row.names(allOfit_scale) <- allOfit_norm$scaf

#
p.mat <- cor.mtest(all_the_wins_scale)
all_the_wins_scale_M<-cor(all_the_wins_scale)

corrplot(all_the_wins_scale_M, method="color",  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE ,tl.cex = 0.75 ,number.cex = 0.6
)
#

allOWin_pca <- PCA(all_the_wins_comp[,-c(1,2,25)],graph = F)

allOWin_pca_labs <- cbind.data.frame(allOWin_pca$ind$coord,all_the_wins_comp[c(1:2)])
allOWin_pca_contribs<- cbind.data.frame("Var"=row.names(allOWin_pca$var$coord),allOWin_pca$var$coord)

#allOfit_pca$var$coord
#allOfit_pca$var$contrib

ggplot() +
  geom_point(data=allOWin_pca_labs,aes(x=Dim.1,y=Dim.2,color=chrom)) + 

  theme_bw() +

#  scale_colour_gradient(low = "skyblue", high = "darkgreen", na.value = NA) +
  labs(color="Size (mb)",size="log(Size)",
       x=paste("PC1 (",round(allOWin_pca$eig[1,2],digits=2),"%)",sep=""),y=paste("PC2 (",round(allOWin_pca$eig[2,2],digits=2),"%)",sep="")) #+
  
  ggrepel::geom_text_repel(
    min.segment.length = 0.5,
    data = allOfit_pca_contribs,
    mapping = aes(x=Dim.1*5,y=Dim.2*5,label = Var),
    position = position_dodge(1.5)
  )

  ####centromeres####
  
 # top_hic <- hic_bypos_connect %>% group_by(chrom) %>% top_n(2,mean_connect)
  
  top_sat <- All_repeats_scafs_sub %>% filter(family == "Satellite")  %>% group_by(chrom) %>% top_n(2,proportion_of_window)  #summarise(mb=mb,highLTRwin = max(proportion_of_window))
  
  bottom_rec <- windowedreC_link %>% group_by(LG) %>% top_n(-2,cmmb)  #summarise(mb=mb,highLTRwin = max(proportion_of_window))
  
  #left_join(top_LTR[,c(1,2,28)] , bottom_rev[,c(3,7,8)],by=c("chr"="LG"))
  
  top_sat<- top_sat[,c(1,2,30)]  
  names(top_sat) <- c("chr","mb","value")
  top_sat$metric <- "topSat"
  top_sat <- 
    top_sat %>% arrange(chr, value) %>%
    group_by(chr) %>% 
    mutate(rank = rank(value))
  
  #bottom_rec <- bottom_rec[,c(3,7,8)]
  names(bottom_rec) <- c("chr","mb","value")
  bottom_rec$metric <- "low cM/mb"
  bottom_rec <- 
    bottom_rec %>% arrange(chr, value) %>%
    group_by(chr) %>% 
    mutate(rank = rank(value))
  
  
 # top_hic <- top_hic[,c(10,5,4)]
 # names(top_hic) <- c("chr","mb","value")
 # top_hic$metric <- "top Hi-C"
 # top_hic <- 
 #   top_hic %>% arrange(chr, value) %>%
 #   group_by(chr) %>% 
 #   mutate(rank = rank(value))
  
  
  
  find_cent <- 
    rbind.data.frame(top_sat,bottom_rec)
  
  find_cent <- left_join(find_cent,sizes,by=c("chr"="chrom"))
  
 # find_cent$Chrom <- factor(find_cent$chr, size_order$chrom)
  find_cent$Chrom <- factor(find_cent$chr, levels=c("1"  ,"1A","2" , "3"  ,"4" ,"4A", "5","Z",  "6",  "7",  "8" , "9" ,"10" ,"11", "12" ,"13", "14", "15" , "17" ,"18", "19", "20" ,"21","22", "23", "24", "25", "26", "27", "28","29"))
  
  
  top_plot_cent <- 
  ggplot(data=find_cent[find_cent$chromType == "macro",]) +
    geom_rect( mapping=aes(xmin=0, xmax=(chrom_length/1000000)+1, ymin=1, ymax=2), fill="white",color="black") + 
   # geom_rect( mapping=aes(xmin=mb, xmax=mb+1, ymin=1, ymax=2, fill=metric)) + 
    geom_point( mapping=aes(x=mb,  y=1.5,  color=metric),size=5) + 
    
      geom_text( mapping=aes(x=mb,  y=1.5,  label=rank)) + 
    
    guides(alpha=FALSE,fill=FALSE,color=FALSE) +
    
    facet_grid(~Chrom,space="free",scales = "free") +  labs(x="Position (mb)") +
    theme(strip.background =element_rect(fill="white"),axis.text.y = element_blank(),axis.ticks.y =  element_blank()) + #xlim(120,125) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())   +
    theme(panel.spacing = unit(0.1, "lines"))
  
  int_plot_cent <- 
    ggplot(data=find_cent[find_cent$chromType == "int",]) +
    geom_rect( mapping=aes(xmin=0, xmax=(chrom_length/1000000)+1, ymin=1, ymax=2), fill="white",color="black") + 
  #  geom_rect( mapping=aes(xmin=mb, xmax=mb+1, ymin=1, ymax=2, fill=metric,alpha=-rank)) + 
    geom_point( mapping=aes(x=mb,  y=1.5,  color=metric),size=5) + 
    
    geom_text( mapping=aes(x=mb,  y=1.5,  label=rank)) + 
    
     guides(alpha=FALSE,fill=FALSE,color=FALSE) +
    facet_grid(.~Chrom,space="free",scales = "free") +  labs(x="Position (mb)") +
    theme(strip.background =element_rect(fill="white"),axis.text.y = element_blank(),axis.ticks.y =  element_blank()) + #xlim(120,125) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())   +
    theme(panel.spacing = unit(0.1, "lines"))
  
  micro_plot_cent <- 
    ggplot(data=find_cent[find_cent$chromType == "micro",]) +
    geom_rect( mapping=aes(xmin=0, xmax=(chrom_length/1000000)+1, ymin=1, ymax=2), fill="white",color="black") + 
 #   geom_rect( mapping=aes(xmin=mb, xmax=mb+1, ymin=1, ymax=2, fill=metric,alpha=-rank)) + 
    geom_point( mapping=aes(x=mb,  y=1.5,  color=metric),size=5) + 
    
    geom_text( mapping=aes(x=mb,  y=1.5,  label=rank)) + 
    
      guides(alpha=FALSE) +
    
    facet_grid(.~Chrom,space="free",scales = "free") +  labs(x="Position (mb)") +
    theme(strip.background =element_rect(fill="white"),axis.text.y = element_blank(),axis.ticks.y =  element_blank()) + #xlim(120,125) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())   +
    theme(panel.spacing = unit(0.1, "lines"))
  
  plot_grid(top_plot,  int_plot, micro_plot, ncol = 3)
  

  
  plot_grid(top_plot_TE,  top_plot_marey, top_plot_cent,ncol = 1,align = 'v',axis='rl')
  plot_grid(but_plot_TE,  but_plot_marey, int_plot_cent,ncol = 1,align = 'v',axis='rl')
  plot_grid(micro_plot_TE,  micro_plot_marey, micro_plot_cent,ncol = 1,align = 'v',axis='rl')
  
  