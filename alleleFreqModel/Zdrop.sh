#!/bin/bash
#SBATCH --partition=rosalind  --time=10-00:00:00  --output=Zdrop.r220315.log
#SBATCH -c 20 --mem=200G 

module load R
module load sendemail/1.56

echo "Z sample"
Rscript --vanilla alleleFreqModelZ_sample.R
sendEmail -f fbeaudry@ur.rochester.edu -t fbeaudry@ur.rochester.edu -u Zsample -m Zsample

echo "Z sim"
Rscript --vanilla alleleFreqModelZ_sim.R
sendEmail -f fbeaudry@ur.rochester.edu -t fbeaudry@ur.rochester.edu -u Zsim -m Zsim

#echo "Z boot"
#Rscript --vanilla alleleFreqModelZ_boot.R