 #!/bin/bash

#pipeline for estimating expected genetic contributions of recent immigrants
#input files: FSJpedgeno_Zsexlinked.ped & coreDemonestlings1990.txt
#Rose Driscoll
#Last updated: 28 September 2021

#generate ped files for contribution of all immigrants or cohorts of immigrants (yearly)
R -f GenContrib_Imm_input.R --vanilla

#run gene dropping
#if you haven't already compiled geneDrop_RD.c, compile now: gcc geneDrop_RD.c -lgsl
##assigning unsexed individuals a sex, with varying sex ratio

for sexRatio in 0 1 0.48
do
    for chrom in A Z 
    do
        python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribAll_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt ImmContribAll_prop_male_${sexRatio}_${chrom} 1000000 0 s ${chrom} ${sexRatio}
    done
done
