#!/bin/bash

#pipeline for estimating expected genetic contributions of 926 individual breeders
#input files: FSJpedgeno_Zsexlinked.ped & IndivDataUSFWS.txt
#Nancy Chen, Rose Driscoll and Felix Beaudry
#Last updated: 18 Aug 2021

#generate ped files for the individuals
R -f GenContrib_Indiv_input.r --vanilla

#run gene dropping sims - specify number of processors to parallelize over in nprocs
#if you haven't already compiled geneDrop_RD.c, compile now: gcc geneDrop_RD.c -lgsl
i=0
nprocs=10 #number of processors to use

for j in IndivContrib_*.ped
do
    python geneDrop_runner_RD.py "$j" working_files/intermediate_files/allABSnestlings.txt "$j".A 1000000 0 s A 1 &
    ((i++))
    [[ $((i%nprocs)) -eq 0 ]] && wait
done
wait

for j in IndivContrib_*.ped
do
    python geneDrop_runner_RD.py "$j" working_files/intermediate_files/allABSnestlings.txt "$j".Z 1000000 0 s Z 1 &
    ((i++))
    [[ $((i%nprocs)) -eq 0 ]] && wait
done
wait