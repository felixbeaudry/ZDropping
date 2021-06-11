#Jun 11 2021
#Rose Driscoll & Felix Beaudry
#script to plot the genetic contributions of male and female immigrants for autosomes and Z

library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(forcats)
library(plyr)
library(fitdistrplus)
library(vcd)
library(foreach)
#setwd('~/Google Drive/Research/Data2/fsj/databases')
#load('all_tables_v2.rdata')

setwd('~/Google Drive/Research/Data2/fsj/ZDropping_all')

#import
#load("simindivFIXmin2obs.rdata")
#load("FSJpedgeno_A.rdata")

# load Rdata file 
load("plotImmGenContrib_tidy_20190529.Rdata")
pedigree_immdata_count <- ungroup(pedigree_immdata_count)

# code to produce all of these tables from raw data is in plotImmGenContrib_tidy_20190418.R

# Nancy's plot theme
plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3), 
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"), 
                    plot.background = element_rect(fill = "white"),  
                    axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), 
                    axis.title = element_text(size=7), plot.title = element_text(size=8), 
                    legend.position="right", legend.text = element_text(size=7),
                    legend.title = element_text(size=8), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(1,"cm"))


####demo_plot####
# Read in and tidy pedigree
ped <- read.table("FSJpedgeno_Zsexlinked.ped", header = FALSE, sep = " ", stringsAsFactors = FALSE)
pedigree <- ped[,1:6]
colnames(pedigree) <- c("Fam", "Indiv", "Dad", "Mom", "Sex", "Pheno")
# Read in immigrant cohort data
indivdata <- read.table("IndivDataUSFWS.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# just need immigrant cohorts
immdata <- dplyr::select(indivdata, Indiv, ImmCohort)


# And join the immigrant cohorts to the pedigree
pedigree_immdata <- left_join(pedigree, immdata, by = "Indiv")

pedigree_immdata_tally <- pedigree_immdata %>% group_by(Sex,ImmCohort) %>% tally()

pedigree_immdata_ratios <- foreach(year= c(1992:2012),.combine=rbind) %do% {
#for(year in c(1992:2012)){
  #  print(year)
  tmp <- pedigree_immdata_tally[pedigree_immdata_tally$ImmCohort == year,]
  res <- chisq.test(c(na.omit(tmp$n[tmp$Sex == 1]),na.omit(tmp$n[tmp$Sex == 2])), p = c(0.5,0.5))
  cbind(year,"p"=res$p.value,"f"=(na.omit(tmp$n[tmp$Sex == 2])/(na.omit(tmp$n[tmp$Sex == 2])+na.omit(tmp$n[tmp$Sex == 1])) ),"t"=(na.omit(tmp$n[tmp$Sex == 2])+na.omit(tmp$n[tmp$Sex == 1])),"u"=na.omit(tmp$n[tmp$Sex == 0]))
}

pedigree_immdata_ratios <- as.data.frame(pedigree_immdata_ratios)
length(pedigree_immdata_ratios$p[pedigree_immdata_ratios$f > 0.5])
length(pedigree_immdata_ratios$p[pedigree_immdata_ratios$p < 0.05])
min(pedigree_immdata_ratios$f)
max(pedigree_immdata_ratios$f)
mean(pedigree_immdata_ratios$f)
min(pedigree_immdata_ratios$t)
max(pedigree_immdata_ratios$t)
mean(pedigree_immdata_ratios$u)

pedigree_immdata_imms_perc <- 
  pedigree_immdata_imms_count %>%
  group_by( ImmCohort) %>% transmute(Sex, percent = num_imms/sum(num_imms))
ggplot(pedigree_immdata_imms_perc,aes(x=ImmCohort,y=percent,fill=as.factor(Sex))) + 
  geom_bar( stat='identity') +
 # facet_grid(does_breed ~ .) + 
  theme_bw()


# Figure out which imms end up breeding
pedigree_immdata_breeder_nonbreeder_imms <- filter(pedigree_immdata, !is.na(ImmCohort)) %>%
  mutate(does_breed = ifelse(Indiv %in% pedigree_immdata$Dad, TRUE, 
                             ifelse(Indiv %in% pedigree_immdata$Mom, TRUE, FALSE)))
# Now count up how many imms that end up breeding or not come in each year
pedigree_immdata_breeder_nonbreeder_imms_count <- group_by(pedigree_immdata_breeder_nonbreeder_imms, ImmCohort, Sex, does_breed) %>%
  dplyr::summarize(num_imms = n()) %>%
  complete(ImmCohort, Sex, does_breed, fill = list(num_imms = 0)) %>%
  filter(!is.na(ImmCohort))
pedigree_immdata_breeder_nonbreeder_imms_perc <- pedigree_immdata_breeder_nonbreeder_imms_count %>%
  group_by( ImmCohort,does_breed) %>% transmute(Sex, percent = num_imms/sum(num_imms))


ggplot(pedigree_immdata_breeder_nonbreeder_imms_perc,aes(x=ImmCohort,y=percent,fill=as.factor(Sex))) + 
  geom_bar( stat='identity') +
  facet_grid(does_breed ~ .) + theme_bw()

pedigree_immdata_breeder_ratios <- foreach(year=c(1991:2012),.combine=rbind) %do% {
  tmp <- pedigree_immdata_breeder_nonbreeder_imms_count[pedigree_immdata_breeder_nonbreeder_imms_count$ImmCohort == year & pedigree_immdata_breeder_nonbreeder_imms_count$does_breed == TRUE,]
  res <- chisq.test(c(tmp$num_imms[tmp$Sex == 1],tmp$num_imms[tmp$Sex == 2]), p = c(0.5,0.5))
  cbind(year,"p"=res$p.value,"f"=(tmp$num_imms[tmp$Sex == 2]/(tmp$num_imms[tmp$Sex == 2]+tmp$num_imms[tmp$Sex == 1]) ),"t"=(tmp$num_imms[tmp$Sex == 2]+tmp$num_imms[tmp$Sex == 1]))
}


pedigree_immdata_breeder_ratios <- as.data.frame(pedigree_immdata_breeder_ratios)
length(pedigree_immdata_breeder_ratios$p[pedigree_immdata_breeder_ratios$f > 0.5])
length(pedigree_immdata_breeder_ratios$p[pedigree_immdata_breeder_ratios$p < 0.05])

min(pedigree_immdata_breeder_ratios$f)
max(pedigree_immdata_breeder_ratios$f)
mean(pedigree_immdata_breeder_ratios$f)

min(pedigree_immdata_breeder_ratios$t)
max(pedigree_immdata_breeder_ratios$t)

ggplot(pedigree_immdata_breeder_ratios,aes(x=year,y=f)) + geom_point()

year_immi_lm <- lm(year~t*f,data=pedigree_immdata_breeder_ratios)
summary(year_immi_lm)


# make males negative
pedigree_immdata_breeder_nonbreeder_imms_count_males_negative <- mutate(pedigree_immdata_breeder_nonbreeder_imms_count, num_imms_males_negative = ifelse(Sex==1, -num_imms, num_imms))

# perhaps a combo plot
# plot sexed birds
(sexed_breed_nonbreed <- ggplot(filter(pedigree_immdata_breeder_nonbreeder_imms_count_males_negative, Sex != 0)) +
    geom_bar(aes(x = ImmCohort, y= num_imms_males_negative, fill = interaction(does_breed, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("lightskyblue1", "cornflowerblue", "mistyrose", "indianred1")) +
    labs(x = "", y = "Number of immigrants"
         #, title = "Immigrants arriving each year that breed by 2013 (dark shades) or do not breed (light shades, gray)"
         ) +
    plottheme +
    scale_y_continuous(breaks = c(-30, -20, -10, 0, 10, 20, 30, 40, 50), labels = c(30, 20, 10, 0, 10, 20, 30, 40, 50)) +
    theme(legend.position = "none", plot.margin=unit(c(0.2,0.1,0,0.15),'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()))
# plot unsexed birds
(unsexed_breed_nonbreed <- ggplot(filter(pedigree_immdata_breeder_nonbreeder_imms_count_males_negative, Sex == 0)) +
    geom_bar(aes(x = ImmCohort, y= num_imms_males_negative, fill = interaction(does_breed, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("gray85")) +
    labs(x = "Year", y = "") +
    plottheme +
    scale_y_continuous(breaks = c(0, 10, 20)) +
    theme(legend.position = "none", plot.margin=unit(c(0.2,0.1,0,0.15),'cm')))

# plot to make a legend
for_legend <- mutate(pedigree_immdata_breeder_nonbreeder_imms_count_males_negative, fill_group = paste(does_breed, Sex, sep = "."))
for_legend$fill_group <- factor(for_legend$fill_group, levels = c("FALSE.2", "TRUE.2", "TRUE.1", "FALSE.1", "FALSE.0"))
(plot_for_legend <- ggplot(for_legend) +
    geom_bar(aes(x = ImmCohort, y= num_imms_males_negative, fill = fill_group), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("mistyrose", "indianred1", "cornflowerblue", "lightskyblue1", "gray85"), labels = c("Female\nnonbreeder", "Female\nbreeder", "Male\nbreeder", "Male\nnonbreeder", "Unsexed\nnonbreeder")) +
    labs(fill = "") +
    plottheme +
    theme(legend.key.size = unit(0.7,"cm")))
sexed_unsexed_breed_nonbreed_legend <- get_legend(plot_for_legend)
# combine
sexed_unsexed_breed_nonbreed__plot <- plot_grid(sexed_breed_nonbreed, unsexed_breed_nonbreed, ncol = 1, nrow = 2, rel_heights = c(0.8, 0.25), align = "v")
#demoplot <- plot_grid(sexed_unsexed_breed_nonbreed__plot, sexed_unsexed_breed_nonbreed_legend, ncol = 2, rel_widths = c(1, 0.2))

# Autosomal contributions

A_imm_gen_contribs <- 
  ggplot(simsumImmA_from1990) + 
  geom_ribbon(aes(x = year, ymin = q1, ymax = q3, fill = as.factor(allele)), alpha=0.3) +
  geom_line(aes(x = year, y = mean, col = as.factor(allele)), size=0.3) + 
  labs(x = "Year", y = "Autosomal expected genetic contrib.") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("cornflowerblue", "indianred1")) +
  scale_fill_manual(values = c("cornflowerblue", "indianred1")) +
  plottheme + 
  theme(legend.position='none',plot.margin=unit(c(0.2,0.1,0,0.15),'cm'))

# Z contributions

Z_imm_gen_contribs <- 
  ggplot(simsumImmZ_from1990) + 
  geom_ribbon(aes(x = year, ymin = q1, ymax = q3, fill = as.factor(allele)), alpha=0.3) +
  geom_line(aes(x = year, y = mean, col = as.factor(allele)), size=0.3) + 
  labs(x = "Year", y = "Z expected genetic contrib.") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("cornflowerblue", "indianred1")) +
  scale_fill_manual(values = c("cornflowerblue", "indianred1")) +
  plottheme + 
  theme(legend.position='none',plot.margin=unit(c(0.2,0.1,0,0.15),'cm'))


####compare A and Z ####
names(simsumImmZ_from1990) <- c("year"  , "allele" ,"mean_Z" ,  "q1_Z"   ,  "q3_Z"    )
names(simsumImmA_from1990) <- c("year"  , "allele" ,"mean_A" ,  "q1_A"   ,  "q3_A"    )

simsumImm_from1990 <- left_join(simsumImmZ_from1990,simsumImmA_from1990,by=c("year"="year","allele"="allele"))
simsumImm_from1990$allele[simsumImm_from1990$allele == 2] <- "M"
simsumImm_from1990$allele[simsumImm_from1990$allele == 3] <- "F"

AZ_immicontrib_plot <-
ggplot(simsumImm_from1990,aes(y=mean_Z,x=mean_A)) +
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  geom_abline(slope = (2/1)*(1/3),  linetype = "dashed", color = "indianred1") +
  geom_abline(slope = (2/1)*(2/3),  linetype = "dashed", color = "cornflowerblue") +
  geom_point(size=2,aes(color=as.factor(allele))) +
  scale_color_manual(values = c( "indianred1","cornflowerblue"))+
  scale_fill_manual(values = c("indianred1","cornflowerblue"))+
  guides(color=FALSE)+
 # geom_smooth(method="lm",color="black") +
 # geom_smooth(method="lm",aes(color=as.factor(allele))) +
  geom_errorbar(aes(ymax =  q3_Z, ymin =  q1_Z,color=as.factor(allele))) +
  geom_errorbarh(aes(xmax =  q3_A, xmin =  q1_A,color=as.factor(allele))) +
  labs(x = "Mean A Immigrant Contrib.", y = "Mean Z Immigrant Contrib.") +plottheme 
 # theme_bw(base_size = 18)+
 # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))


ggplot(simsumImm_from1990,aes(y=mean_Z,x=mean_A,color=year,shape=as.factor(allele))) +
  geom_abline(slope = 1,  linetype = "dashed", color = "gray") +
  geom_abline(slope = (2/1)*(1/3),  linetype = "dashed", color = "indianred1") +
  geom_abline(slope = (2/1)*(2/3),  linetype = "dashed", color = "cornflowerblue") +
  geom_point(size=2) +
  geom_errorbar(aes(ymax =  q3_Z, ymin =  q1_Z)) +
  geom_errorbarh(aes(xmax =  q3_A, xmin =  q1_A)) +
  labs(x = "Mean A", y = "Mean Z",title="Cumul. Contrib. of Immigrants") +
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))

AZ_mod <- lm(mean_A ~ 0 + mean_Z , data=simsumImm_from1990)  
summary(AZ_mod)  

AZ_mod <- lm(mean_Z ~ 0 + mean_A , data=simsumImm_from1990)  
summary(AZ_mod)  

simsumImm_from1990$ZA <-  0.5 * (simsumImm_from1990$mean_Z/simsumImm_from1990$mean_A) 

ggplot(simsumImm_from1990,aes(x=year,y=ZA,color=as.factor(allele))) + geom_point() + geom_line() + stat_smooth() +
  theme_bw(base_size = 22) + scale_color_manual(values=c(   '#00ADEF',   '#00A550'  )) +labs(x="Year",y="Z/A E(H)",color="Sex")



####Full Plotting ####

plot_grid(sexed_unsexed_breed_nonbreed__plot, A_imm_gen_contribs,AZ_immicontrib_plot, Z_imm_gen_contribs,ncol = 2,labels = c("A","C","B","D"), label_size = 12,align = 'hv',axis='tbrl')


ggplot(simsumImm_from1990,aes(y=mean_Z/mean_A,x=year)) +
  annotate('segment', x=1993,xend=2013,y=1,yend=1,alpha=0.1)+
  
  annotate('segment', x=1993,xend=2013,y=(4/3),yend=(4/3),alpha=0.5,linetype="dashed")+
  annotate('text',x=1991,y=4/3,label='4/3',size=5)+

  annotate('segment', x=1993,xend=2013,y=(2/3),yend=(2/3),alpha=0.5,linetype="dashed")+
  annotate('text',x=1991,y=2/3,label='2/3',size=5)+
  
  
    geom_line(aes(color=as.factor(allele)))+
  geom_ribbon(aes(x = year, ymin = q1_Z/mean_A, ymax = q3_Z/mean_A, fill = as.factor(allele)), alpha=0.3) +
  
  geom_ribbon(aes(x = year, ymin = q1_A/mean_A, ymax = q3_A/mean_A, fill = as.factor(allele)), alpha=0.3) +
  
  labs(x = "Year", y = "Z/A Immigrant Contributions",color="Sex",fill="Sex") +
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.margin = unit(c(15,5.5,5.5,5.5), "pt"))


