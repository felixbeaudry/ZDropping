# Felix Beaudry & Nancy Chen
# format input file for SNPeff

# rename chromosomes with the full scaffold name
# load in file with snp position information 
snp.pos <- read.table("snp.pos.txt", header = TRUE, stringsAsFactors=FALSE)

# import vcf
library(data.table)
library(tidyverse)
realvcf <- fread("geno.anon.vcf",skip = "#")

#remove individuals
realvcf9 <- realvcf[,c(1:9)]
realvcf9_names <- names(realvcf9) #save names for later

#join mapped markers
fauxvcf <- left_join(realvcf9,snp.pos,by=c("ID"="SNP"))

#remove old columns
fauxvcf_reorder <- fauxvcf[,c("scaffold","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]

# edit unmapped SNPs
fauxvcf_reorder[fauxvcf_reorder$scaffold == '*', 'scaffold'] <- 0

#sort by chromosome and position
fauxvcf_sorted <- fauxvcf_reorder[order(fauxvcf_reorder$scaffold,fauxvcf_reorder$POS),]

#give back old names
names(fauxvcf_sorted) <- realvcf9_names

write.table(fauxvcf_sorted, file = "geno.anon.faux.vcf", 
            append = F, quote = FALSE, sep = "\t", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, qmethod = c("escape", "double"))
