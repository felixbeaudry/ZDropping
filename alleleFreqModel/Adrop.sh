#!/bin/bash
#SBATCH --partition=rosalind  --time=10-00:00:00  --output=Adrop.r220128.log
#SBATCH -c 20 --mem=200G  

module load r/4.1.1/b1
module load sendemail/1.56


echo "A sample"
Rscript --vanilla alleleFreqModelA_sample.R
sendEmail -f fbeaudry@ur.rochester.edu -t fbeaudry@ur.rochester.edu -u Adropsample -m Adropsample


echo "A sim"
Rscript --vanilla alleleFreqModelA_sim.R
#echo "A boot"
#Rscript --vanilla alleleFreqModelA_boot.R

sendEmail -f fbeaudry@ur.rochester.edu -t fbeaudry@ur.rochester.edu -u Adropsim -m Adropsim
