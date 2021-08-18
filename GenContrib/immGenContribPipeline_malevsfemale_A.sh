#!/bin/bash

#pipeline for estimating expected genetic contributions of recent immigrants
#input files: FSJpedgeno_Zsexlinked.ped & coreDemonestlings1990.txt
#Rose Driscoll
#Last updated: 3 April 2019

#generate ped files for contribution of all immigrants or cohorts of immigrants (yearly)
#R -f immPed_malevsfemale.r --vanilla

#run gene dropping
#if you haven't already compiled geneDrop_final.c, compile now: gcc geneDrop_final.c -lgsl
python geneDrop_runner_RD.py ImmContribAll_malevsfemale.ped allABSnestlings.txt ImmContribAll_malevsfemale_A 1000000 0 s A
python geneDrop_runner_RD.py ImmContribYearly_malevsfemale.ped allABSnestlings.txt ImmContribYearly_malevsfemale_A 1000000 0 s A
