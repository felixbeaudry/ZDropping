genome_version = FSJ2

#regenerate size file (with the scaffolded assembly)
python juicer/misc/generate_site_positions.py DpnII ${genome_version} ${genome_version}.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${genome_version}_DpnII.txt > genome_files/${genome_version}.scaffSizes.txt

##NUCMER to Zebra Finch##
#SBATCH --mem=100G
 
module load mummer
#gunzip GCA_009859065.2_bTaeGut2.pri.v2_genomic.fna.gz

nucmer -l 1000 -p mummer/${genome_version}_TaeGut \
	${genome_version}.fasta \
	GCA_009859065.2_bTaeGut2.pri.v2_genomic.fna

delta-filter -i 50 -l 1000 -q \
	-r ${genome_version}_TaeGut.delta \
	> ${genome_version}_TaeGut.filti50l1000.delta
show-coords -c -d -l -r \
	 -T ${genome_version}_TaeGut.filti50l1000.delta \
	 > ${genome_version}_TaeGut.filti50l1000.delta.coords

##REPEATS##
#SBATCH -c 21 --mem=50G

module load repeatmodeler genometools/1.5.9 mafft/7.313 cdhit/4.6.8

BuildDatabase -name "${genome_version}" ${genome_version}.fa

RepeatModeler -database ${genome_version} -pa 7 -LTRStruct \
	-genometools_dir genometools/1.5.9/bin/  \
	-cdhit_dir cdhit/4.6.8/ \
	-mafft_dir mafft/7.313/bin/ \
	-ltr_retriever_dir LTR_retriever/ \
	-ninja_dir NINJA/NINJA 

#stick Valentina's library out repeatmodeler's -- they don't have same width (fold)
awk '{print $1}' ${genome_version}-families.fa >${genome_version}.repeatLib.fa
fold -w 50 genome_files/peona_bird_library_25Oct2020.lib >>${genome_version}.repeatLib.fa

mkdir ${genome_version}.mask
RepeatMasker -s -lib genome_files/${genome_version}.repeatLib.fa -dir ${genome_version}.mask ${genome_version}.hifiasm.scaff.fasta -pa 7 -xsmall

##JUICER WITH REPEAT MASKED ASSEMBLY##

module load gcc/9.1.0 jre/1.8.0_111 bwa/0.7.4 juicer/1.5.6 python3

bwa index ${genome_version}.fa.masked

#this is a repeat of L2 and L3, but in a different place
python juicer/misc/generate_site_positions.py DpnII ${genome_version} ${genome_version}.fa.masked
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${genome_version}_DpnII.txt > ${genome_version}.chrom.sizes

##input directly into command line, sbatch cannot handle this!
for library in 402 403
do
	juicer/scripts/juicer_local.sh -d juicer/DTG-HiC-${library} -D juicer -z juicer/references/${genome_version}.fa.masked -p juicer/restriction_sites/FSJ3.chrom.sizes -s juicer/restriction_sites/${genome_version}_DpnII.txt -s DpnII -q ${server_name} -l ${server_name}  -L 180:00:00 -Q 180:00:00 -C 180000000 -g ${genome_version} 
done

sh juicer/scripts/mega_local.sh -g ${genome_version} -s DpnII -D juicer -q ${server_name} -l ${server_name}  -L 180:00:00 -Q 180:00:00

java -Xms512m -Xmx2048m -jar juicer_tools_1.22.01.jar pre  -q 30 -r 1000000 mega/aligned/merged_nodups.txt mega/aligned/inter_Q30_res.1mb.hic ${genome_version}.chrom.sizes

awk '$2 > 1000000 {print $1}' restriction_sites/${genome_version}.chrom.sizes >longscafs.list


while read chrom1
do
	while read chrom2 
	do

		java -Xms512m -Xmx2048m -jar juicer_tools_1.22.01.jar dump \
			observed NONE inter_Q30_res.1mb.hic \
			${chrom1} ${chrom2} BP 1000000 \
			${chrom1}_${chrom2}.1mb.raw.hic.txt

		awk -v chrom1=$chrom1 -v chrom2=$chrom2 '{print chrom1"\t"chrom2"\t"$1"\t"$2"\t"$3}' \
		${chrom1}_${chrom2}.1mb.raw.hic.txt >>${genome_version}.hic.1mb.raw.txt

	done < juicer/longscafs.list
done < juicer/longscafs.list

##fixed distance between SNPs from crimap

while read LG 
do
	crimap chr${LG}.par fixed -n 36 > chr${LG}.fixed.txt

	sed -e '1,/Sex_averaged/ d' chr${LG}.fixed.txt | sed '/denote/q' | sed '/denote/d' | tail -n +2 | paste -d " "  - -  | awk -v chr=chr${chr} '{print chr"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"}' >>crimap.fixed.txt

done < LG.list



##REMAP CHIP MARKERS on repeatmasked##

module load perl bwa
bwa index ${genome_version}.fasta.gz
bwa mem ${genome_version}.fasta.gz beadchipSeq.fq >${genome_version}.beadChip.sam 

perl linkerScripts/seqLocationSNP.pl ${genome_version}.beadChip.sam ${genome_version}.beadChip.sam.raw.pos 

sed -e 's/caffold_//g' ${genome_version}.beadChip.sam.raw.pos >${genome_version}.beadChip.sam.pos

#ldhat methods by Tram tn337@cornell.edu, then reformatted + pruned using linkerScripts/pruneLDhat.R


##CG AND CPG ISLANDS##

module load perl samtools python3

#mkdir scafs
while read chrom
do
	echo "CpGs of ${chrom}"

	samtools faidx ${genome_version}.fa ${chrom} >scafs/${chrom}.fa
	perl linkerScripts/CpGcluster.pl scafs/${chrom}.fa 50 1E-5 scafs/${chrom}.fasta.cpg

done < ${genome_version}.scafs

awk '$1 ~ "CGI" {print "chrom\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' scafs/Chromosome1.fasta.cpg >scafs/all.cpg

while read chrom
do 
	
	awk -v chrom=$chrom '$1 !~ "CGI" {print chrom"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' scafs/${chrom}.fasta.cpg >>scafs/all.cpg
	
done < /scratch/fbeaudry/V3.scafs

python3 linkerScripts/gapper.py -i ${genome_version}.fa -w 1000 >genome_files/${genome_version}.gapsINwin1mb.txt

