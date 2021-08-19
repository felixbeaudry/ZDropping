 #!/bin/bash

#pipeline for estimating expected genetic contributions of recent immigrants
#input files: FSJpedgeno_Zsexlinked.ped & coreDemonestlings1990.txt
#Rose Driscoll
#Last updated: 3 April 2019

#generate ped files for contribution of all immigrants or cohorts of immigrants (yearly)
R -f GenContrib_Imm_input.R --vanilla

#run gene dropping
#if you haven't already compiled geneDrop_final.c, compile now: gcc geneDrop_final.c -lgsl
python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribAll_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribAll_malevsfemale_A 1000000 0 s A
python geneDrop_runner_RD.py ImmContribYearly_malevsfemale.ped allABSnestlings.txt ImmContribYearly_malevsfemale_A 1000000 0 s A

python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt ImmContribAll_malevsfemale_Z 1000000 0 s Z
python geneDrop_runner_RD.py ImmContribYearly_malevsfemale.ped allABSnestlings.txt ImmContribYearly_malevsfemale_Z 1000000 0 s Z


##assigning unsexed individuals a sex, with varying sex ratio

#run gene dropping for Z & autosomes with proportion male 1 (all unsexed individuals assigned male)
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_1_Z 1000000 0 s Z 1
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_1_A 1000000 0 s A 1

#run gene dropping for Z & autosomes with proportion male 0 (all unsexed individuals assigned female)
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_0_Z 1000000 0 s Z 0
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_0_A 1000000 0 s A 0

#run gene dropping for Z & autosomes with proportion male 1 (all unsexed individuals assigned male)
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_0.48_Z 1000000 0 s Z 0.48
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_0.48_A 1000000 0 s A 0.48
