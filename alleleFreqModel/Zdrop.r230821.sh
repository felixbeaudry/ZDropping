#!/bin/bash
#SBATCH --partition=rosalind  --time=10-00:00:00  --output=Zdrop.r230821.log
#SBATCH 

module load r/4.0.5/b1

Rscript --vanilla alleleFreqModelZ_sample.R
Rscript --vanilla alleleFreqModelZ_sim.R
Rscript --vanilla alleleFreqModelZ_boot.R

module load sendemail/1.56
sendEmail -f fbeaudry@ur.rochester.edu -t fbeaudry@ur.rochester.edu -u Zdrop -m Zdrop



module load sendemail/1.56
sendEmail -f fbeaudry@ur.rochester.edu -t fbeaudry@ur.rochester.edu -u Adrop -m Adrop
