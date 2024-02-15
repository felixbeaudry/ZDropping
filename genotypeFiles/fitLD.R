####ld####

ld.bychr <- fread('pedigree.bychr.ld')[,-c(1,2,4,5)] %>% left_join(marey_zerod[,c(1,13,12)],by=c("SNP_A"="SNPname")) %>% left_join(marey_zerod[,c(1,13,12)],by=c("SNP_B"="SNPname")) %>% filter(homolog.x == homolog.y)

ld.bychr$mb_distance <- abs(ld.bychr$mb.x - ld.bychr$mb.y)

ld.bychr$R2 <- as.numeric(ld.bychr$R2)

ld.bychr$chrType[ld.bychr$homolog.x != "Z"] <- "A"
ld.bychr$chrType[ld.bychr$homolog.x == "Z"] <- "Z"

fit_A <- lm(na.omit(ld.bychr$R2[ ld.bychr$chrType == "A"]) ~ I(1/na.omit(ld.bychr$mb_distance[ ld.bychr$chrType == "A"])))
summary(fit_A)

ld_hyperbola_A <- function(x){
  y = fit_A$coefficients[2]/x + fit_A$coefficients[1]
  print(y)
}

ld_hyperbola_A(1) - ld_hyperbola_A(5) 

fit_Z <- lm(na.omit(ld.bychr$R2[ ld.bychr$chrType == "Z"]) ~ I(1/na.omit(ld.bychr$mb_distance[ ld.bychr$chrType == "Z"])))

ld_hyperbola_Z <- function(x){
  y = fit_Z$coefficients[2]/x + fit_Z$coefficients[1]
}


pdf(paste("fig_S2.pdf",sep=''),width=5.5,height=4)

ggplot(ld.bychr ,aes(x=mb_distance,y=R2,color=as.factor(chrType))) +  
  geom_point(size=0.5,alpha=0.01) + 
  theme_bw()  + ylim(0,1)+ xlim(0,0.1)+ 
  geom_function(fun = ld_hyperbola_A, colour = "red") +
  geom_function(fun = ld_hyperbola_Z, colour = "blue") + 
  scale_color_manual(values=c("red","blue")) +
  labs(x="Pairwise Distance (mb)",
       y=expression(paste("Linkage Disequilibrium (",R^2,")")),color="Chrom. Type") +
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))

dev.off()