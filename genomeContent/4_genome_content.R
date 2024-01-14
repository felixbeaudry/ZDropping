#genome_content by Felix Beaudry, April 9th 2022

setwd('~/Documents/Github/genomeContent/')
source("genomeContent_functions.R")

####start####

window_size = 1000000

genomeName="FSJ2"

primaryAssembly <- fread(paste0('genome_files/',genomeName,'.scaffSizes.txt'))

#summarize sizes
sizes_primary <- assemblySizes(primaryAssembly)

#make one file with scaff name, length in bp and bird homolog (homolog_up and _down are just factors in ascending/descending order for plotting)
size_homology <- homologyMap(conversionFile=paste0('genome_files/',genomeName,'.convert.txt'),size=sizes_primary) #Chicken-based homology

#split each scaffold into wins (of window_size)
df_map <- windowedMap(sizes=sizes_primary,window_size = window_size) 


#makes dataframe of GC content and number of Ns/gaps & plots GC content
bybase <- calcGapsNBases(bybaseFile=paste0('genome_files/',genomeName,'.gapsINwin',window_size/1000000,'mb.txt'),
                         homology_conversion=size_homology,
                         gap_window= window_size)

#makes dataframe of CpG islands and plots it
cpg <- makeCpG(cpgFile=paste0('genome_files/',genomeName,'.cpg.txt'),homology_conversion=size_homology) ##cpg islands##

#makes a dataframe of repeat content per window by repeat family
repeats_windowed <- makeRepeatWindows(genomeName=genomeName,
                                      df_map=df_map,
                                      homology_conversion=size_homology,
                                      window_size_TE=window_size)

#makes dataframe of HiC links per chromosome
hic_byChrom <- hic_summarizer(genomeName,window_size,size_homology)

#makes dataframe of coding/non-coding-expressed bp per window
coding_windows <- coding_windower(annotation_winData=paste0('genome_files/',genomeName,'.annotation_win.Rdata'), 
                                  size_homology,window_size)


####Recombination####

genmap <- makeGeneticMap(genmapFile='genome_files/crimap.fixed.txt')

marey <- makeMareyMap(chipFile=paste0('genome_files/',genomeName,'.beadChip.sam.pos'),genmap=genmap,size_homology=size_homology)

rec_windows_cm <- get_rec_windows(CMMB = marey %>% dplyr::select(homolog,SNPname,cm_flip,mb),wSize = window_size/1000000) #columns = "chr"    "marker" "cm"     "mb"    


ldhat_processed <- ldhatFit(ldhatFile=paste0('genome_files/',genomeName,'.rho.Rdata'),size_homology=size_homology)
rec_windows_rho <- get_rec_windows(CMMB = ldhat_processed %>% dplyr::select(homolog,SNPID,Cum_rho_adj,mb)) # cols: "chr", "marker", "cm", "mb"

#compare linkage and rho
names(rec_windows_cm)[3] <- "link"
names(rec_windows_rho)[3] <- "rho"

windowedreC_join <- left_join(rec_windows_cm,rec_windows_rho)

windowedreC_corr <- compareRHO2MB(windowedreC_join=windowedreC_join,size_homology=size_homology)

#sex & linkage
sex_marey <- makeSexMarey(sexlinkagefile='genome_files/crimap.fixed.sex.txt',marey=marey)

rec_windows_male <- get_rec_windows(CMMB = sex_marey %>% filter(variable == "male_adj") %>% dplyr::select(lg,marker,cm_flip,mb) ) # cols: "lg", "marker", "cm", "mb"
rec_windows_female <- get_rec_windows(CMMB = sex_marey %>% filter(variable == "female_adj") %>% dplyr::select(lg,marker,cm_flip,mb) ) # cols: "lg", "marker", "cm", "mb"

rec_SDI <- makeLinkageSDI(rec_windows_male=rec_windows_male,rec_windows_female=rec_windows_female,marey=marey)
#when male is large, value is negative

#fit hyperbola to LD data
fit_ld(marey,ld_file = 'genome_files/ld.txt.gz')

####Sequence Content####
#makes dataframe of HiC links per chromosome
hic_byChrom <- hic_summarizer(genomeName,window_size,size_homology)

#makes a dataframe of repeat content per window by repeat family
repeats_windowed <- makeRepeatWindows(genomeName=genomeName,
                                      df_map=df_map,
                                      homology_conversion=size_homology,
                                      window_size_TE=window_size)

#makes dataframe of coding/non-coding-expressed bp per window
coding_windows <- coding_windower(annotation_winData=paste0('genome_files/',genomeName,'.annotation_win.Rdata'), 
                                  size_homology,
                                  window_size,
                                  genomeName=genomeName)

#makes dataframe of GC content and number of Ns/gaps & plots GC content
bybase <- calcGapsNBases(bybaseFile=paste0('genome_files/',genomeName,'.gapsINwin',window_size/1000000,'mb.txt'),
                         homology_conversion=size_homology,
                         gap_window= window_size)

#makes dataframe of CpG islands and plots it
cpg <- makeCpG(cpgFile=paste0('genome_files/',genomeName,'.cpg.txt'),homology_conversion=size_homology) ##cpg islands##


####combine variables into one dataframe####

repeats_grouped <- 
  repeats_windowed %>% filter(TEfam %in% c("DNA_transposon","house_RNAs","LINE","Low_complexity","LTR","Non-LTR","Satellite","Simple_repeat")) %>% 
  group_by(homolog,TEfam) %>% summarize(win_mean=mean(proportion_of_window)) %>% dcast(homolog~TEfam)

hic_byChrom <-  hic_byChrom %>% group_by(chr_1) %>% summarise(hic=mean(meanHiC))

GCContent_byChrom <- bybase %>% group_by(homolog) %>% summarise(meanGC = mean(GC))

rec_SDI_byChrom <- rec_SDI %>% group_by(homolog) %>% summarise(meanRECSDI = mean(SDI,na.rm=T))


allOfit <-  
  
  left_join(size_homology %>% filter(!is.na(homolog)), GCContent_byChrom ) %>%
  left_join(cpg %>% dplyr::select(homolog,CpG)) %>%
  left_join(repeats_grouped) %>%
  left_join(rec_SDI_byChrom) %>%
  left_join(hic_byChrom,by=c("homolog"="chr_1")) %>%
#  left_join(coding_windows %>% group_by(homolog,label) %>% summarise(ave_bp=mean(fracCoding)) %>% dcast(formula=homolog~label)) %>%
  left_join(marey %>% group_by(homolog) %>% summarise(map_length=max(cm_flip))) %>%
  left_join(sex_marey %>% group_by(homolog,variable_name) %>% summarise(map_length=max(cm_flip)) %>% dcast(formula=homolog~variable_name)) 

allOfit <- allOfit %>% filter(bp > 2000000)

row.names(allOfit) <- allOfit$homolog

save(allOfit,file=paste0('genome_files/',genomeName,'.genomeContent.Rdata'))




