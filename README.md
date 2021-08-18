# ZDropping
Sex-biased demography, including sex-biased survival or migration, can impact allele frequency changes across the genome. In particular, we expect different patterns of genetic variation on autosomes and sex chromosomes due to sex-specific differences in life histories, as well as differences in effective population size, transmission modes, and the strength and mode of selection. Here, we present a set of scripts to directly characterize the relative roles of sex-biased demography and inheritance in shaping genome-wide allele frequency trajectories. These scripts are associated with Driscoll et al. (biorXiv).

Scripts are seperated into directories roughly corresponding to the sections of the paper and, replicating the order of the methods in the above paper. 

## Base Dropping
The first step for the majority of subsequent analyses is the gene dropping simulation. The scripts for the basic gene dropping step is in the directory `\base_Dropping`. The step is performed by running the python script `geneDrop_runner_RD.py` which is a wrapper for geneDrop_RD.c (compiled as a.out using `gcc geneDrop_RD.c -lgsl`). To run `geneDrop_runner_RD.py [ped file] [cohort file] [prefix for output files] [number of drops] [filter ungenotyped individuals? Y/N] [output type] [A, Z, or X] [proportion of unsexed individuals to assign as male]`

## Expected Genetic Contributions
The next step in our procedure was to analyze expected genetic contributions for individuals and for immmigrants as a group. The scripts to support these analyses are in the directory `\GenContrib`. The scripts calculate and plot expected genetics contributions given the pedigree. This includes scripts to calculate expected contributions for autosomal or Z_linked loci, as well as scripts to calculate expected contributions for individuals or for all immigrant individuals. 

### Individual expected genetic contributions
We first calculate expected genetic contributions for each individual in ABS using scripts in the subdirectory `\individual`. To do this, we run `GenContrib_Indiv.sh` which will first call `GenContrib_Indiv_input.r` to make the input pedigrees for each individual and then pass these pedigree to `geneDrop_runner_RD.py` first using an autosomal model of inheritance, then using a Z model of inheritance. 
With expected genetic contributions for each individual, we first illustrate principals of expected genetic contributions using one chosen couple from the pedigree. `GenContrib_plot_fig1_couple.R` plots figure 1. Next, we plot all individual expected genetic contributions using `GenContrib_Indiv_plot_fig2.R`.

### Immigrant expected genetic contributions
We next calculate expected genetic contributions for all immigrants to ABS using scripts in the subdirectory `\immigrant`. The pipeline in `GenContrib_Imm.sh` will first make the appropriate input files using `GenContrib_Imm_input.R`, then calling `geneDrop_runner_RD.py` first using an autosomal model of inheritance, then using a Z model of inheritance. Finally, GenContrib_Imm.sh runs through `geneDrop_runner_RD.py` while varing the sex ratio assigned to unsexed individuals.
With expected genetic contributions for all immigrants, we can plot figure 3. We first run `plotImmGenContrib_tidy_20190418.R` to calculate the expected genetic contributions of immigrant male and females then `GenContrib_Imm_plot_fig3_runDE.R` to compare the expected genetic contributions of immigrants to their cohort size. FInally we merge all datasets into one plot using `GenContrib_Imm_plot_fig3_merge.R`


## Signals Of Selection
Next, we look for signals of selection by comparing expected contributions to observed change in frequency. The directory `\SignalsOfSelection` holds scripts to track actual genotypes down the pedigree for Z-linked loci. The script `analyzeSelOutput_Z.R` will perform the analysis of the results from 

## Allele Frequency Change Model
The final directory `\alleleFreqModel` holds scripts to partition change in allele frequencies between years between demographic groups by sex for autosomal and Z-linked loci. This analysis is broken down into  four steps: 1 - sampling, 2 - simulating error, 3 - bootstrap and 4 - plot. This is repeated for the Z and autosomes seperately, given different inheritance patterns, but both are plotted together. These scripts assume a flexible sex ratio; we also ran a model assuming a fixed sex ratio but instead partitioned allele frequency change between years to each sex. 

### Sample
Step 1 - run `alleleFreqModelA_sample.R` and `alleleFreqModelZ_sample.R` to calculate variance in allele frequencies between years. These scripts will also setup the files necessary for the next step.

### Simulate error
Step 2 - run `alleleFreqModelA_sim.R` and `alleleFreqModelZ_sim.R` for autosomal and Z loci respectively to calculate error.

### Bootstrap

Step 3 - run `alleleFreqModelA_boot.R` and `alleleFreqModelZ_boot.R` to calculate confidence intervals 

### Plot

Step 4 - `alleleFreqModelAZ_plot.R` will combine Z and autosomal analyses into one set of plots for easy comparison.

### Fixed Sex Ratio
If you are interested in the fixed sex ratio model, these scripts are in the directory `\fixedSexRatio_model`, and follow the same format as steps 1, 2 and 4 above; we do not run bootstraps on the fixed sex ratio model.