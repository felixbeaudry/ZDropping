# genomeContent

## Table of Contents
1. [Dependencies](#dependencies)
2. [Step by Step walkthrough](##step-by-step)

## Introduction
`genomeContent` can be thought of as a suite of recommendations to annotate various aspects of a new genome. While mostly relying on pre-existing software and packages, `genomeContent` has some recommended parameters to smooth over this somewhat disparate literature. Please cite all the appropriate dependencies below. To note, this suite was developped for the annotation of the genome of the Florida Scrub-Jay, so makes some assumptions biased towards the assembly of bird genome such as widespread karyotype homology leading to pre-existing chromosome numbering.

## Dependencies
### Software
- [bwa](http://bio-bwa.sourceforge.net/) v.0.7.4
- [juicer](https://github.com/aidenlab/juicer) v.1.5.6
	- [gcc](https://gcc.gnu.org/) 9.1.0
	- [jre](https://java.com/en/download/) 1.8.0_111
- [Chromonomer](http://catchenlab.life.illinois.edu/chromonomer/) v1.13
- [nucmer](http://mummer.sourceforge.net/) v.4.0.0
- [repeatmodeler](https://www.repeatmasker.org/RepeatModeler/) v.2.0.1
- [repeatmasker](https://www.repeatmasker.org/RepeatModeler/) v.4.1.0
	- [genometools](http://genometools.org/) v1.5.9
	- [mafft](https://mafft.cbrc.jp/alignment/software/) v7.313
	- [cdhit](http://weizhong-lab.ucsd.edu/cd-hit/) 4.6.8
	- [LTR_retriever](https://github.com/oushujun/LTR_retriever) v2.9.0
	- [ninja](https://github.com/TravisWheelerLab/NINJA) v0.98
- [busco](https://busco.ezlab.org/) v.5.2.2
- [BRAKER](https://github.com/Gaius-Augustus/BRAKER) v.2.1.6
	- [diamond](http://github.com/bbuchfink/diamond) v.0.9.24
	- [prothint](https://github.com/gatech-genemark/) v.2.6.0
	- [augustus](https://github.com/Gaius-Augustus/Augustus) v3.4.0
	- [genemark](http://exon.gatech.edu/GeneMark/license_download.cgi) v4.68
	- [cbdfasta](https://github.com/gpertea/cdbfasta) (sync'd on 2021-01-18)
- [picard](https://broadinstitute.github.io/picard/) v.2.12.0
- [samtools](http://www.htslib.org/) v.1.11
- [stringtie](https://ccb.jhu.edu/software/stringtie/)
- [star](https://github.com/alexdobin/STAR) v.2.7.3
- [bedtools](bedtools.readthedocs.org) v.2.30
- [snpeff](http://pcingola.github.io/SnpEff/se_running/) v4.3

## Step by Step
First, run `1_Assembly.sh`. We further scaffold the Dovetail-scaffolded assembly using a preliminary linkage map. We first map markers to the new assembly using `bwa`. We then make a dummy agp file using the chromosome lengths (which I generate using the `generate_site_positions.py` script from `juicer` - but there's likely a much easier way - followed by a little in-house R script called `linkerScripts/make_agp.R`). `Chromonomer` then makes a new .agp informed by the linkage map, and `agp2assembly.py` by `juicer` to put the assembly together into a fasta file. There's then a few lines to add scaffolds missed by chromonomer and to change the names of the scaffolds, if you'd like.

Next, run `2_gene_annotation.sh` to annotate the genome. First, we use `nucmer` et al. to find homology to Zebra finch. Then we use `repeatmodeler` to make de novo repeat annotation; cat that with the annotation from [Peona et al. 2020](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13252) to make a repeat library and run through `repeatmasker`. Then, it'll run `busco` for the whole genome to get of a sense of how complete the genome is (we'll also be using the whole genome busco file as a sort of scaffolding for the other annotation). Then, we make a `BRAKER2` run (the protein version). To seed the protein annotation, we stick together the vertebrate orthodb database and all the chicken genes on unipro. Run `busco` on that. Next `BRAKER1` uses RNAseq data, so first we align with `star`, then `sort` reads using `picard` and `samtools`. Run `busco` on that too. Finally, we annotate using `stringtie` (that step is technically missing from the script). Run `busco` on that. Because I couldn't find a script that I liked to combine these things, I wrote `geneCombo.R` which will combine all these outputs into a windowed analysis of the genes of the genome as well as output an annotation that has all the best versions of the genes identified in each of the above independent annotations.

Finally, run `3_SNP_annotation.sh` to create a VCF file and annotate the SNPs using SNPeff.

In the directory, `/genome_files` you'll find examples for some of the files you need. 

## Funding
The development of `genomeContent` was supported by a National Science Foundation (NSF) Postdoctoral Research Fellowship in Biology, grant number 2109639. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NSF.

