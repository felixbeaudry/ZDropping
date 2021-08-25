# ZDropping
Sex-biased demography, including sex-biased survival or migration, can impact allele frequency changes across the genome. In particular, we expect different patterns of genetic variation on autosomes and sex chromosomes due to sex-specific differences in life histories, as well as differences in effective population size, transmission modes, and the strength and mode of selection. Here, we present a set of scripts to directly characterize the relative roles of sex-biased demography and inheritance in shaping genome-wide allele frequency trajectories. These scripts are associated with Driscoll et al. (biorXiv). This work also build on [Chen et al 2019](https://www.pnas.org/content/116/6/2158) with associated scripts and data at [here](http://dx.doi.org/10.6084/m9.figshare.7044368)

Scripts are seperated into directories roughly corresponding to the sections of the paper and, replicating the order of the methods in the above paper. 

## Expected Genetic Contributions
For this section, analyses are nested across several scripts, with slight differences for the two main subsections: Individual and Immigrant expected genetic contributions. The general workflow is in the bash script `GenContrib_*.sh`, and proceed through making the intermediate files using `GenContrib_*_input.R`, which requires two input files: a raw pedigree file (eg `working_files\pedigree.txt`) and raw information about cohort file (eg `working_files\IndivDataUSFWS.txt`). The output files are the input for the next step.

The next step is running the python script `geneDrop_runner_RD.py` which is a wrapper for geneDrop_RD.c (compiled as a.out using `gcc geneDrop_RD.c -lgsl`; Requires [GSL](https://www.gnu.org/software/gsl/doc/html/) which you can install with `brew install gsl`). To run `geneDrop_runner_RD.py` set the following parameters in this order [processed ped file; string] [processed cohort file; string] [prefix for output files] [number of drops; int] [filter ungenotyped individuals? 0/1] [output type; s] [A, Z, or X; string] [proportion of unsexed individuals to assign as male; int]; for example `python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribAll_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribAll_malevsfemale_A 1000000 0 s A 1`.

Finally, we plot results using the `GenContrib_*_plot_*.R` scripts. These plotting scripts are labelled by their associated figures (`_fig*`).



### Signals Of Selection
Next, we look for signals of selection by comparing expected contributions to observed change in frequency. The directory `\SignalsOfSelection` holds scripts to track actual genotypes down the pedigree for Z-linked loci. The script `analyzeSelOutput_Z.R` will perform the analysis of the results from 

## Allele Frequency Change Model
The final directory `\alleleFreqModel` holds scripts to partition change in allele frequencies between years between demographic groups by sex for autosomal and Z-linked loci. This analysis is broken down into  four steps: 1 - sampling, 2 - simulating error, 3 - bootstrap and 4 - plot. This is repeated for the Z and autosomes seperately, given different inheritance patterns, but both are plotted together. These scripts assume a flexible sex ratio; we also ran a model assuming a fixed sex ratio but instead partitioned allele frequency change between years to each sex. 

### Step 1 - Sample
run `alleleFreqModelA_sample.R` and `alleleFreqModelZ_sample.R` to calculate variance in allele frequencies between years. These scripts will also setup the files necessary for the next step.

### Step 2 - Simulate error
 run `alleleFreqModelA_sim.R` and `alleleFreqModelZ_sim.R` for autosomal and Z loci respectively to calculate error.

### Step 3 - Bootstrap

run `alleleFreqModelA_boot.R` and `alleleFreqModelZ_boot.R` to calculate confidence intervals 

### Step 4 - Plot

`alleleFreqModelAZ_plot.R` will combine Z and autosomal analyses into one set of plots for easy comparison.

### Fixed Sex Ratio model
If you are interested in the fixed sex ratio model, these scripts are in the directory `\fixedSexRatio_model`, and follow the same format as steps 1, 2 and 4 above; we do not run bootstraps on the fixed sex ratio model.