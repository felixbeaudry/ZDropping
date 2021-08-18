#!/bin/bash

#pipeline for estimating expected genetic contributions of recent immigrants
#input files: FSJpedgeno_Zsexlinked.ped & IndivDataUSFWS.txt
#Nancy Chen and Rose Driscoll
#Last updated: 17 Apr 2018

#generate ped files for contribution of all immigrants or cohorts of immigrants (yearly)
#R -f indivPed_Z.r --vanilla

#run gene dropping sims - specify number of processors to parallelize over in nprocs
#if you haven't already compiled geneDrop_final.c, compile now: gcc geneDrop_final.c -lgsl
i=0
nprocs=10 #number of processors to use

for j in IndivContrib_*.ped
do
    python geneDrop_runner_RD.py "$j" allABSnestlings.txt "$j".Z 1000000 0 s Z 1 &
    ((i++))
    [[ $((i%nprocs)) -eq 0 ]] && wait
done
wait
