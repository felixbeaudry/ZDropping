# script to make manhattan plots of the selection tests data
# Rose Driscoll
# Last updated: 04 October 2021


## Setup

library(ggplot2)
library(dplyr)
library(cowplot)

# load in data tables
load("working_files/intermediate_files/SignalsOfSelection_Results.rdata")
# will also need autosomal results for multiple comparisons corrections
load("working_files/SignalsOfSelection_autosomalResults.rdata")

# new Z snp order for plotting & tables
Z_SNP_order <- read.table("working_files/Z_SNP_order.txt", header = TRUE, stringsAsFactors = FALSE)

# add new Z snp order to tables
pval_1999to2013 <- left_join(pval_1999to2013, Z_SNP_order, by = c("snp" = "ped_order"))
pval_change <- left_join(pval_change, Z_SNP_order, by = c("snp" = "ped_order"))

# plot theme
plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3), 
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"), 
                    plot.background = element_rect(fill = "white"),  
                    axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), 
                    axis.title = element_text(size=7), plot.title = element_text(size=8), 
                    legend.position="right", legend.text = element_text(size=7),
                    legend.title = element_text(size=8), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(1,"cm"))


## Adjust p-values and assess significance

#Apply pseudocounts to get rid of 0s - we will recalculate pvals, adding 1 to the # of sims </> the obs
pval_1999to2013 <- mutate(pval_1999to2013, pval_new = ((pval*1000000)+1)/1000000)
pval_change <- mutate(pval_change, pval_new = ((pval*1000000)+1)/1000000)
#Also apply pseudocounts to autosomal data
pval_1999to2013_autosomes <- mutate(pval_1999to2013_autosomes, pval_new = ((pval_median*1000000)+1)/1000000)
pval_change_autosomes <- mutate(pval_change_autosomes, pval_new = ((pval_change_median*1000000)+1)/1000000)

# Using FDR method, adjust for all snps in genome (268 Z snps + 10731 autosomal snps)

# adjust 1999-2013 comparisons
pval_1999to2013$padjAll <- p.adjust(c(pval_1999to2013$pval_new, pval_1999to2013_autosomes$pval_new), method = "fdr", n = 268+10731)[1:268]
# Apply cutoff
pval_1999to2013 <- mutate(pval_1999to2013, padjAll_sig = ifelse(padjAll < 0.25, 2, 1))
filter(pval_1999to2013, padjAll_sig==2)

# adjust year-by-year comparisons
collect_padj <- NULL
pval_change_w_padj <- NULL
for (i in 1999:2012) {
  collect_padj <- p.adjust(c(pval_change[pval_change$year == i, 'pval'], pval_change_autosomes[pval_change_autosomes$year == i, 'pval_new']), method = "fdr", n = 268+10731)[1:268]
  pval_change_thisyear <- filter(pval_change, year ==i)
  pval_change_thisyear$padjAll <- collect_padj
  pval_change_w_padj <- rbind(pval_change_w_padj, pval_change_thisyear)
  }
pval_change_w_padj <- mutate(pval_change_w_padj, padjAll_sig = ifelse(padjAll < 0.25, 2, 1))
filter(pval_change_w_padj, padjAll_sig==2)

# Check whether the autosomal significance results have changed

# 1999-2013
pval_1999to2013_autosomes$padjAll <- p.adjust(c(pval_1999to2013$pval_new, pval_1999to2013_autosomes$pval_new), method = "fdr", n = 268+10731)[269:10999]
filter(pval_1999to2013_autosomes, adjpval_median1 < 0.25, padjAll >= 0.25)
# year-by-year
collect_padj <- NULL
pval_change_w_padj_auto <- NULL
for (i in 1999:2012) {
  collect_padj <- p.adjust(c(pval_change[pval_change$year == i, 'pval'], pval_change_autosomes[pval_change_autosomes$year == i, 'pval_new']), method = "fdr", n = 268+10731)[269:10999]
  pval_change_thisyear <- filter(pval_change_autosomes, year ==i)
  pval_change_thisyear$padjAll <- collect_padj
  pval_change_w_padj_auto <- rbind(pval_change_w_padj_auto, pval_change_thisyear)
}
filter(pval_change_w_padj_auto, adjpval_median < 0.25, padjAll >= 0.25)


## Plot

# plot with significant 1999-2013 snps after adjusting for all snps in the genome highlighted in orange
ggplot(data = pval_1999to2013) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  scale_x_continuous(name='Position on Z chromosome') +
  scale_y_continuous(name='-log10(p-value)') +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  labs(title = "Significant SNPs 1999-2013")
#ggsave("manhattanplots_all.pdf", width = 4, height = 3, units = "in")

# 1999-2000
(manhattan_1999_2000 <- ggplot(data = filter(pval_change_w_padj, year == 1999)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = 'Significant SNPs by year\n1999-2000') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2000-2001
(manhattan_2000_2001 <- ggplot(data = filter(pval_change_w_padj, year == 2000)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '\n2000-2001') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2001-2002
(manhattan_2001_2002 <- ggplot(data = filter(pval_change_w_padj, year == 2001)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2001-2002') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2002-2003
(manhattan_2002_2003 <- ggplot(data = filter(pval_change_w_padj, year == 2002)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2002-2003') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2003-2004
(manhattan_2003_2004 <- ggplot(data = filter(pval_change_w_padj, year == 2003)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2003-2004') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2004-2005
(manhattan_2004_2005 <- ggplot(data = filter(pval_change_w_padj, year == 2004)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2004-2005') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2005-2006
(manhattan_2005_2006 <- ggplot(data = filter(pval_change_w_padj, year == 2005)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2005-2006') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2006-2007
(manhattan_2006_2007 <- ggplot(data = filter(pval_change_w_padj, year == 2006)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2006-2007') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2007-2008
(manhattan_2007_2008 <- ggplot(data = filter(pval_change_w_padj, year == 2007)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2007-2008') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2008-2009
(manhattan_2008_2009 <- ggplot(data = filter(pval_change_w_padj, year == 2008)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2008-2009') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2009-2010
(manhattan_2009_2010 <- ggplot(data = filter(pval_change_w_padj, year == 2009)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2009-2010') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2010-2011
(manhattan_2010_2011 <- ggplot(data = filter(pval_change_w_padj, year == 2010)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2010-2011') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2011-2012
(manhattan_2011_2012 <- ggplot(data = filter(pval_change_w_padj, year == 2011)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2011-2012') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

# 2012-2013
(manhattan_2012_2013 <- ggplot(data = filter(pval_change_w_padj, year == 2012)) + 
  geom_point(aes(x = SNP_order, y = -log10(pval_new), color=as.factor(padjAll_sig)), alpha = 0.75, size = 0.3) +
  plottheme + 
  labs(x='Position on Z chromosome', y = '-log10(p-value)', title = '2012-2013') +
  scale_y_continuous(limits = c(0,6)) +
  scale_colour_manual(values=c("black","#fc8d59")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()))

plot_grid(manhattan_1999_2000, manhattan_2000_2001, manhattan_2001_2002, manhattan_2002_2003, manhattan_2003_2004, manhattan_2004_2005, manhattan_2005_2006, manhattan_2006_2007, manhattan_2007_2008, manhattan_2008_2009, manhattan_2009_2010, manhattan_2010_2011, manhattan_2011_2012, manhattan_2012_2013, nrow = 7, ncol = 2, rel_heights = c(1.1,1,1,1,1,1,1))
#ggsave("manhattanplots_yearly_6-5x8-5.pdf", width = 6.5, height = 8.5, units = "in")
