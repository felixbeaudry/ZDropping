# ZDropping
Sex-biased demography, including sex-biased survival or migration, can impact allele frequency changes across the genome. In particular, we expect different patterns of genetic variation on autosomes and sex chromosomes due to sex-specific differences in life histories, as well as differences in effective population size, transmission modes, and the strength and mode of selection. Here, we present a set of scripts to directly characterize the relative roles of sex-biased demography and inheritance in shaping genome-wide allele frequency trajectories. These scripts are associated with Driscoll et al. (biorXiv). We seperate our scripts into three directory, replicating the order of the methods in the above paper. 

The first directory 'expectedContributions' holds the scripts to calculate and plot expected genetics contributions given the pedigree. This includes scripts to calculate expected contributions for autosomal or Z_linked loci, as well as scripts to calculate expected contributions for individuals or for immigrant individuals as a whole.

The second directory 'SignalsOfSelection' holds scripts to track actual genotypes down the pedigree for Z-linked loci.

The third directory 'alleleFreqModel' holds scripts to partition change in allele frequencies between years between demographic groups by sex for autosomal and Z-linked loci.