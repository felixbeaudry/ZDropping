#!/bin/bash

#pipeline for using gene-dropping to test for short-term selection
#input files: FSJpedgeno_Zsexlinked.ped, FSJpedgeno_Zpseudoautosomal.ped & IndivDataUSFWS.txt
#Rose Driscoll
#Last updated: 03 October 2021

#generate cohort file for genedropping simulations
R -f SignalsOfSelection_coreCohort.R --vanilla

#if you haven't already compiled geneDrop_RD.c, compile now: gcc geneDrop_RD.c -lgsl

python geneDrop_runner_RD.py working_files/FSJpedgeno_Zsexlinked.ped working_files/coreDemoNestlings.txt working_files/intermediate_files/seltest_Zsexlinked 1000000 1 l Z 1
python geneDrop_runner_RD.py working_files/FSJpedgeno_Zpseudoautosomal.ped working_files/coreDemoNestlings.txt working_files/intermediate_files/seltest_Zpseudoautosomal 1000000 1 l A 1

