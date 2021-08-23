#!/bin/bash
#SBATCH --partition=rosalind  --time=10-00:00:00  --output=Adrop.r230821.log
#SBATCH 

module load r/4.0.5/b1

Rscript --vanilla alleleFreqModelA_sample.R
Rscript --vanilla alleleFreqModelA_sim.R
Rscript --vanilla alleleFreqModelA_boot.R