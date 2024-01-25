library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
FSJ3.chrom.sizes <- fread(args[1]) 

agp <- 
  cbind.data.frame(
    FSJ3.chrom.sizes$V1,
    1,
    FSJ3.chrom.sizes$V2,
    1,
    "W",
    FSJ3.chrom.sizes$V1,
    1,
    FSJ3.chrom.sizes$V2,
    "?"
  )

write.table(agp, file = "FSJ3.hifiasm.agp", 
            append = F, quote = FALSE, sep = "\t", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))