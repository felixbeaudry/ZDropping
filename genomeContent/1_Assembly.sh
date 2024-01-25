genome_version = FSJ2

##SECONDARY ASSEMBLY WITH CHROMONOMER##
module load bwa perl R samtools python

#index
bwa index ${genome_version}.fa

#align chip markers to primary assembly
bwa mem ${genome_version} \
	genome_files/beadchipSeq.fq >${genome_version}.beadChip.sam

#make list of scaffold lengths
python juicer/misc/generate_site_positions.py DpnII ${genome_version} \
	${genome_version}.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${genome_version}_DpnII.txt > ${genome_version}.chrom.sizes

#makes an agp from sizes file, only designed for primary assembly
Rscript --vanilla linkerScripts/make_agp.R ${genome_version}.chrom.sizes

chromonomer-1.13/chromonomer -p genome_files/crimap.fixed.txt  \
	-o ./ \
	-s ${genome_version}.beadChip.sam  \
	-a ${genome_version}.agp \
	--data_version ${genome_version}

#chromonomer outputs an agp-type file, next two steps makes new fasta from that agp
#https://github.com/phasegenomics/juicebox_scripts/blob/master/juicebox_scripts
python juicebox_scripts-master/juicebox_scripts/agp2assembly.py \
	${genome_version}.agp \
	${genome_version}.agp.assembly

python juicebox_scripts-master/juicebox_scripts/juicebox_assembly_converter.py \
	-f ${genome_version}.fa \
	-a ${genome_version}.agp.assembly \
	-p ${genome_version}.scaff.fa -v 

#the agp from chromonomer has no info about contigs not scaffolded, these awk steps add those back into assembly fasta
awk '{print $1}' ${genome_version}.chrom.sizes | \
	sort >${genome_version}.chrom
awk '{print $6}' ${genome_version}.agp |\
	 sort | uniq | tail -n +2 >scaffolded_contigs.txt

#comm compares two files
comm -23 ${genome_version}.chrom \
	scaffolded_contigs.txt \
	>unscaffolded_contigs.txt

cp ${genome_version}.scaff.fa.fasta ${genome_version}.fa
echo -e "\n\n" >>${genome_version}.fa
#stick on remaining scaffolds
samtools faidx --mark-strand no ${genome_version}.fa -n 80 \
	-r unscaffolded_contigs.txt \
	 >>${genome_version}.fa
