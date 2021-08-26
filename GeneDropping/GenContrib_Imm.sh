 #!/bin/bash

#pipeline for estimating expected genetic contributions of recent immigrants
#input files: FSJpedgeno_Zsexlinked.ped & coreDemonestlings1990.txt
#Rose Driscoll
#Last updated: 3 April 2019

#generate ped files for contribution of all immigrants or cohorts of immigrants (yearly)
R -f GenContrib_Imm_input.R --vanilla

#run gene dropping
#if you haven't already compiled geneDrop_RD.c, compile now: gcc geneDrop_RD.c -lgsl
python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribAll_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribAll_malevsfemale_A 1000000 0 s A
python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribYearly_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribYearly_malevsfemale_A 1000000 0 s A

python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribAll_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribAll_malevsfemale_Z 1000000 0 s Z
python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribYearly_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribYearly_malevsfemale_Z 1000000 0 s Z


##assigning unsexed individuals a sex, with varying sex ratio

for sexRatio in 0 1 0.48
do
    for chrom in A Z 
    do
        python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_${sexRatio}_${chrom} 1000000 0 s ${chrom} ${sexRatio}
    done
done
