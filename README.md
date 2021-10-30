# ZDropping
Sex-biased demography, including sex-biased survival or migration, can impact allele frequency changes across the genome. In particular, we expect different patterns of genetic variation on autosomes and sex chromosomes due to sex-specific differences in life histories, as well as differences in effective population size, transmission modes, and the strength and mode of selection. Here, we present a set of scripts to directly characterize the relative roles of sex-biased demography and inheritance in shaping genome-wide allele frequency trajectories. These scripts are associated with [Driscoll et al. (biorXiv)](https://www.biorxiv.org/content/10.1101/2021.10.28.466320v1). This work also builds on [Chen et al 2019](https://www.pnas.org/content/116/6/2158), with associated scripts and data [here](http://dx.doi.org/10.6084/m9.figshare.7044368).

Scripts are separated into one directory for those sections of the methods which employ gene dropping (expected genetic contributions of individuals and immigrants, signals of selection) and one section for the allele frequency variance models. An additional preliminary script and directory of genotype files are used to generate input files for downstream analyses.

## Make genotype input files
Raw genotypes are stored in a zipped format to save space. First, unzip the genotype file in the `\genotypeFiles` directory using:

`gzip -d geno.anon.ped.gz`

To generate input files for gene dropping and the allele frequency models, run `makeGenotypeInput.R`.

## Gene Dropping

The `\GeneDropping` directory holds scripts and data files relevant to expected genetic contributions and signals of selection, analyses which employ gene dropping.

### Basic gene dropping simulation
We implement gene dropping in C++ and supply options for autosomal, ZW, or XY transmission patterns. The gene dropping script is run using a python wrapper as described below.

#### Gene dropping script; how to compile
The gene dropping script is `geneDrop_RD.c`. Compiling this script requires [GSL](https://www.gnu.org/software/gsl/doc/html/) which you can install with `brew install gsl`. Compile as follows to produce the file `a.out`:

`gcc geneDrop_RD.c -lgsl`

#### Python wrapper script
The python wrapper script for gene dropping is `geneDrop_runner_RD.py`. The options for this script are as follows: 

[processed ped file; string] [processed cohort file; string] [prefix for output files] [number of drops; int] [filter ungenotyped individuals? 0/1] [output type; s] [A, Z, or X; string] [proportion of unsexed individuals to assign as male; int]

An example of how to run this script would be: `python geneDrop_runner_RD.py working_files/intermediate_files/ImmContribAll_malevsfemale.ped working_files/intermediate_files/allABSnestlings.txt working_files/intermediate_files/ImmContribAll_malevsfemale_A 1000000 0 s A 1`

Each analysis that uses gene dropping (both expected genetic contributions and signals of selection) includes bash scripts with the correct parameters to run `geneDrop_runner_RD.py` for that analysis.

### Expected genetic contributions of individuals
We use gene dropping to simulate the expected genetic contributions of individuals. This analysis is broken down into two steps: 1 - gene dropping, and 2 - plotting. 

#### Step 1 - Gene dropping
Use `GenContrib_Indiv.sh` to run the gene dropping simulations for a set of 926 breeders. This script also runs `GenContrib_Indiv_input.R` in order to generate the pedigrees that are the input for the gene dropping simulations.

#### Step 2 - Plotting
`GenContrib_Indiv_plot_fig1_pair.R` and `GenContrib_Indiv_plot_fig2.R` plot the individual expected genetic contributions for a selected pair of breeders and the full set of 926 breeders respectively.

### Expected genetic contributions of immigrants
We use gene dropping to simulate the combined expected genetic contributions of immigrants into the population. This analysis is broken down into two steps: 1 - gene dropping, and 2 - plotting. 

#### Step 1 - Gene dropping
Use `GenContrib_Imm.sh` to run the gene dropping simulations for immigrants. This script also runs `GenContrib_Imm_input.R` in order to generate the pedigrees that are the input for the gene dropping simulations.

#### Step 2 - Plotting
`GenContrib_Imm_plot_fig3.R` plots the expected genetic contributions of immigrants.

### Signals of selection
Next, we look for signals of selection by comparing expected contributions to observed change in frequency. This analysis is broken down into three steps: 1 - gene dropping, 2 - performing tests of selection, and 3 - plotting. 

#### Step 1 - Gene dropping
Use `SignalsOfSelection_Sel.sh` to run the gene dropping simulations with real genotypes to simulate the neutral behavior of alleles. This script also runs `SignalsOfSelection_coreCohort.R` to generate cohort data for the gene dropping simulations.

#### Step 2 - Perform tests of selection
Run `SignalsOfSelection_analyzeSelOutput_Z.R` to test for selection based on the gene dropping simulations.

#### Step 3 - Plot
Run `SignalsOfSelection_plot_figS5_figS6.R` to correct selection results for multiple comparisons genome-wide and create manhattanplots to display selection results.

## Allele Frequency Change Model
The final directory, `\alleleFreqModel`, holds scripts to partition change in allele frequencies between years between demographic groups by sex for autosomal and Z-linked loci. This analysis is broken down into four steps: 1 - sampling, 2 - simulating error, 3 - bootstrapping, and 4 - plotting. This is repeated for the Z and autosomes seperately, given different inheritance patterns, but both are plotted together. These scripts assume a flexible sex ratio; we also ran a model assuming a fixed sex ratio but instead partitioned allele frequency change between years  relative to the total number of individuals of each sex. 

#### Step 1 - Sample
Run `alleleFreqModelA_sample.R` and `alleleFreqModelZ_sample.R` to calculate variance in allele frequencies between years. These scripts will also set up the files necessary for the next step.

#### Step 2 - Simulate error
Run `alleleFreqModelA_sim.R` and `alleleFreqModelZ_sim.R` for autosomal and Z loci respectively to calculate error.

#### Step 3 - Bootstrap
Run `alleleFreqModelA_boot.R` and `alleleFreqModelZ_boot.R` to calculate confidence intervals 

#### Step 4 - Plot
`alleleFreqModelAZ_plot.R` will combine Z and autosomal analyses into one set of plots for easy comparison.

### Fixed sex ratio model
Scripts for the fixed sex ratio model follow the same format as steps 1, 2 and 4 above. We do not run bootstraps on the fixed sex ratio model.

#### Step 1 - Sample
Run `alleleFreqModelA_sexRatio_sample.R` and `alleleFreqModelZ_sexRatio_sample.R` to calculate variance in allele frequencies between years. These scripts will also set up the files necessary for the next step.

#### Step 2 - Simulate error
Run `alleleFreqModelA_sexRatio_sim.R` and `alleleFreqModelA_sexRatio_sim.R` for autosomal and Z loci respectively to calculate error.

#### Step 4 - Plot
`alleleFreqModelAZ_sexRatio_plot.R` will combine Z and autosomal analyses into one set of plots for easy comparison.
