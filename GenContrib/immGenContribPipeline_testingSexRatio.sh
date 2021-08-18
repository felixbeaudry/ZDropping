#!/bin/bash

# testing my new genedrop runner that assigns unsexed individuals a sex

#run gene dropping for Z & autosomes with proportion male 1 (all unsexed individuals assigned male)
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_1_Z 1000000 0 s Z 1
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_1_A 1000000 0 s A 1

#run gene dropping for Z & autosomes with proportion male 0 (all unsexed individuals assigned female)
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_0_Z 1000000 0 s Z 0
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_0_A 1000000 0 s A 0

#run gene dropping for Z & autosomes with proportion male 1 (all unsexed individuals assigned male)
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_0.48_Z 1000000 0 s Z 0.48
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt test_ImmContribAll_prop_male_0.48_A 1000000 0 s A 0.48
