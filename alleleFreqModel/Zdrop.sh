#!/bin/bash
#SBATCH --partition=rosalind  --time=10-00:00:00  --output=Zdrop.r220128.log
#SBATCH -c 20 --mem=200G 

module load r/4.1.1/b1

echo "Z sample"
Rscript --vanilla alleleFreqModelZ_sample.R
echo "Z sim"
Rscript --vanilla alleleFreqModelZ_sim.R
#echo "Z boot"
#Rscript --vanilla alleleFreqModelZ_boot.R

module load sendemail/1.56
sendEmail -f fbeaudry@ur.rochester.edu -t fbeaudry@ur.rochester.edu -u Zdrop -m Zdrop
