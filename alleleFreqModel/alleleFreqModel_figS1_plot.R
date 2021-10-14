## Plotting the number of survivors, immigrants, and births each year that are included in the allele freq model
## Rose Driscoll
## last updated 19 May 2021


## Setup

library(dplyr)
library(ggplot2)
library(cowplot)

#plot theme
plottheme <- theme( axis.line.x = element_line(colour="black",size=0.3), axis.line.y = element_line(colour="black",size=0.3),
                    axis.ticks = element_line(colour = "black",size=0.2),
                    axis.text = element_text(colour="black"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"),
                    plot.background = element_rect(fill = "white"),
                    axis.text.x = element_text(size=6), axis.text.y = element_text(size=6),
                    axis.title = element_text(size=7), plot.title = element_text(size=8),
                    legend.position="right", legend.text = element_text(size=7),
                    legend.title = element_text(size=8), legend.key = element_rect(colour=NA,fill=NA), legend.key.size=unit(0.25,"cm"))

#indivlist
#load("working_files/simindivFIXmin2obs.rdata")
#ped<-read.table('working_files/FSJpedgeno_Zsexlinked.ped',header=FALSE,sep=' ',stringsAsFactors=FALSE)
#pedinfo <- ped[,1:5]
#colnames(pedinfo) <- c("Family", "USFWS", "Dad", "Mom", "Sex")
#indivlist <- merge(simindivFIXmin2obs[,1:6],pedinfo[,c(2,5)],by='USFWS')
#indivlist <- indivlist[order(indivlist$Year),]
load(file='working_files/intermediate_files/indivlistgeno_A.rdata')
indivlist <- indivlistgeno_A[,c(1:8)]


## Number of genotyped & not genotyped indivs in each category in each year
counts_to_plot <- group_by(indivlist, Year, Category, Sex, Genotyped) %>%
  dplyr::summarize(num_birds = n()) %>%
  # truncate to start from 1991 as 1990 is just "founder"
  filter(Year>1990) %>%
  # make males negative & make is_genotyped (true/false) column
  mutate(num_birds_males_negative = ifelse(Sex==1, -num_birds, num_birds), is_genotyped = ifelse(Genotyped=="Y", TRUE, FALSE))

# Make survivor, immigrant, and birth tables
survivor_counts_to_plot <- filter(counts_to_plot, Category == "survivor")
immigrant_counts_to_plot <- filter(counts_to_plot, Category == "immigrant")
birth_counts_to_plot <- filter(counts_to_plot, Category == "nestling")



## Plot

# plot to make a legend
for_legend <- mutate(counts_to_plot, fill_group = paste(is_genotyped, Sex, sep = "."))
for_legend$fill_group <- factor(for_legend$fill_group, levels = c("FALSE.2", "TRUE.2", "TRUE.1", "FALSE.1", "FALSE.0"))
(plot_for_legend <- ggplot(for_legend) +
    geom_bar(aes(x = Year, y= num_birds_males_negative, fill = fill_group), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("mistyrose", "indianred1", "cornflowerblue", "lightskyblue1", "gray85"), labels = c("Ungenotyped\nfemale", "Genotyped\nfemale", "Genotyped\nmale", "Ungenotyped\nmale", "Ungenotyped\nunsexed\nindividual")) +
    labs(fill = "") +
    plottheme +
    theme(legend.key.size = unit(0.7,"cm")))
counts_plot_legend <- get_legend(plot_for_legend)

# Survivors

# plot sexed birds
(sexed_survivors <- ggplot(filter(survivor_counts_to_plot, Sex != 0)) +
    geom_bar(aes(x = Year, y= num_birds_males_negative, fill = interaction(is_genotyped, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("lightskyblue1", "cornflowerblue", "mistyrose", "indianred1")) +
    labs(x = "", y = "Number of survivors") +
    plottheme +
    scale_y_continuous(breaks = c(-120,-80,-40, 0, 40, 80, 120), 
                       labels = c(120, 80, 40, 0, 40, 80, 120)) +
    scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010), limits = c(1990,2014)) +
    theme(legend.position = "none", plot.margin=unit(c(0.5,0.1,0,0.15),'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), panel.grid.major = element_line(color = "gray95")))
# plot unsexed birds
(unsexed_survivors <- ggplot(filter(survivor_counts_to_plot, Sex == 0)) +
    geom_bar(aes(x = Year, y= num_birds_males_negative, fill = interaction(is_genotyped, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("gray85")) +
    labs(x = "", y = "") +
    plottheme +
    scale_y_continuous(breaks = c(0, 20), limits = c(0,20)) +
    scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010), limits = c(1990,2014)) +
    theme(legend.position = "none", plot.margin=unit(c(0.1,0.1,0.5,0.15),'cm'), panel.grid.major = element_line(color = "gray95")))
# combine
(survivor_plot <- plot_grid(sexed_survivors, unsexed_survivors, ncol = 1, nrow = 2, rel_heights = c(1, 0.3), align = "v", labels = c("A", "")))

# Immigrants

# plot sexed birds
(sexed_immigrants <- ggplot(filter(immigrant_counts_to_plot, Sex != 0)) +
    geom_bar(aes(x = Year, y= num_birds_males_negative, fill = interaction(is_genotyped, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("lightskyblue1", "cornflowerblue", "mistyrose", "indianred1")) +
    labs(x = "", y = "Number of immigrants        ") +
    plottheme +
    scale_y_continuous(breaks = c(-20, 0, 20),labels = c(20, 0, 20), limits = c(-20,20)) +
    scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010), limits = c(1990,2014)) +
    theme(legend.position = "none", plot.margin=unit(c(0.6,0.1,0,0.15),'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), panel.grid.major = element_line(color = "gray95")))
# plot unsexed birds
(unsexed_immigrants <- ggplot(filter(immigrant_counts_to_plot, Sex == 0)) +
    geom_bar(aes(x = Year, y= num_birds_males_negative, fill = interaction(is_genotyped, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("gray85")) +
    labs(x = "", y = "") +
    plottheme +
    scale_y_continuous(breaks = c(0, 5), limits = c(0,5)) +
    scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010), limits = c(1990,2014)) +
    theme(legend.position = "none", plot.margin=unit(c(0.1,0.1,0.5,0.15),'cm'), panel.grid.major = element_line(color = "gray95")))
# combine
(immigrant_plot <- plot_grid(sexed_immigrants, unsexed_immigrants, ncol = 1, nrow = 2, rel_heights = c(1, 0.6), align = "v", labels = c("B", "")))

# Births

# plot sexed birds
(sexed_births <- ggplot(filter(birth_counts_to_plot, Sex != 0)) +
    geom_bar(aes(x = Year, y= num_birds_males_negative, fill = interaction(is_genotyped, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("lightskyblue1", "cornflowerblue", "mistyrose", "indianred1")) +
    labs(x = "", y = "Number of births") +
    plottheme +
    scale_y_continuous(breaks = c(-120, -80, -40, 0, 40, 80, 120),
    labels = c(120, 80, 40, 0, 40, 80, 120), limits = c(-120,120)) +
    scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010), limits = c(1990,2014)) +
    theme(legend.position = "none", plot.margin=unit(c(0.5,0.1,0,0.15),'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), panel.grid.major = element_line(color = "gray95")))
# plot unsexed birds
(unsexed_births <- ggplot(filter(birth_counts_to_plot, Sex == 0)) +
    geom_bar(aes(x = Year, y= num_birds_males_negative, fill = interaction(is_genotyped, Sex)), position = "stack", stat = "identity", alpha = 0.75) +
    scale_fill_manual(values = c("gray85")) +
    labs(x = "Year", y = "") +
    plottheme +
    scale_y_continuous(breaks = c(0, 40, 80, 120), limits = c(0,120)) +
    scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010), limits = c(1990,2014)) +
    theme(legend.position = "none", plot.margin=unit(c(0.1,0.1,0,0.15),'cm'), panel.grid.major = element_line(color = "gray95")))
# combine
(birth_plot <- plot_grid(sexed_births, unsexed_births, ncol = 1, nrow = 2, rel_heights = c(1, 0.5), align = "v", labels = c("C", "")))



## Combine plots

(SIB_plots <- plot_grid(survivor_plot, immigrant_plot, birth_plot, nrow = 3, rel_heights = c(0.3, 0.2, 0.4)))
plot_grid(SIB_plots, counts_plot_legend, ncol = 2, rel_widths = c(1,0.25))
#ggsave("alleleFreq_model_figS1.pdf", width = 6.5, height = 9, units = "in")

