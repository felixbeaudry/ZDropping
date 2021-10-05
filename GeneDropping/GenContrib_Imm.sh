 #!/bin/bash

#pipeline for estimating expected genetic contributions of recent immigrants
#input files: pedigree.txt & IndivDataUSFWS.txt
#Rose Driscoll
#Last updated: 04 October 2021

#generate ped files for contribution of all immigrants or cohorts of immigrants (yearly)
R -f GenContrib_Imm_input.R --vanilla

#run gene dropping
#if you haven't already compiled geneDrop_RD.c, compile now: gcc geneDrop_RD.c -lgsl

## all immigrants
##assigning unsexed individuals a sex, with a few different sex ratios
for sexRatio in 0 1 0.48
do
    for chrom in A Z 
    do
        python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribAll_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribAll_prop_male_${sexRatio}_${chrom} 1000000 0 s ${chrom} ${sexRatio}
    done
done

## immigrant cohorts (yearly)
for chrom in A Z 
    do
        python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribYearly_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribYearly_malevsfemale_${chrom} 1000000 0 s ${chrom} 1
    done
