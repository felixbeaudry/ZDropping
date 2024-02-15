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

#genome BUSCO, to set base-line

module load busco/5.2.2
busco -m genome \
    -i ${genome_version}.mask/${genome_version}.fasta.masked \
    -o ${genome_version}.hardmasked.busco -l aves_odb10 -f -c 20

##PROTEIN##
wget https://v100.orthodb.org/download/odb10_vertebrata_fasta.tar.gz
tar xvf odb10_vertebrata_fasta.tar.gz
cat vertebrate/Rawdata/* > vert_proteins.fasta
fold -w 60 vert_proteins.fasta > all_proteins.fa

#https://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ano+AND+organism%3A%22Gallus+gallus+%28Chicken%29+%5B9031%5D%22+AND+proteome%3Aup000000539

awk '{print $1}' uniprot_chicken.fasta | \
    sed $'s/[^[:alnum:]\t]//g' | \
    sed ':a;N;$!ba;s/tr/\>/g' | \
    sed ':a;N;$!ba;s/sp/\>/g' >uniprot_Gallus_gallus.fa

fold -w 60 uniprot_Gallus_gallus.fa >>all_proteins.fa

sed $'s/://g' all_proteins.fa | \
    sed $'s/>/>gene/g' >all_proteins.nospe.fa


#rm -r braker/GeneMark-E*

module unload perl/5.24.0/b1
module load braker/2.1.6 diamond/0.9.24/b2 prothint/2.6.0 genemark/b1 augustus/3.4.0 perl/5.32.0/b1

#ID runs, in case you run many times (which you very likely will)

perl braker/2.1.6/scripts/braker.pl --softmasking --cores=21 --gff3 \
    --epmode --species ${genome_version} --useexisting \
    --genome=${genome_version}.fa.masked \
    --prot_seq=all_proteins.nospe.fa \
    --workingdir=BRAKER2/${run} \
    --DIAMOND_PATH=diamond/0.9.24/b2/bin/ \
    --PROTHINT_PATH=prothint/2.6.0/bin/ \
    --GENEMARK_PATH=genemark/4.68 \
    --AUGUSTUS_SCRIPTS_PATH=augustus/3.4.0/scripts \
    --AUGUSTUS_BIN_PATH=augustus/3.4.0/bin \
    --CDBTOOLS_PATH=cdbfasta/ \
    --AUGUSTUS_CONFIG_PATH=augustus/Augustus/config/ 

module unload python3/3.8.3

module load busco/5.2.2 bedtools

awk '$3 ~ "gene" {print $1"\t"$2"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' BRAKER2/${run}/braker.gtf > BRAKER2/${run}/braker.${run}.gtf

bedtools getfasta -fi ${genome_version}.mask/${genome_version}.fa.masked \
    -bed BRAKER2/${run}/braker.${run}.gtf \
    -nameOnly -fo BRAKER2/${run}/braker.${run}.gtf.fa

busco -m transcriptome -i BRAKER2/${run}/braker.${run}.gtf.fa \
    -o ${genome_version}.BRAKER2_busco.tsv -l aves_odb10 -f -c 21


##RNA##

#SBATCH -c 21 --mem=50G

module load rnastar/2.7.3a picard/2.12.0 samtools

STAR --runMode genomeGenerate \
    --genomeDir ${genome_version}.mask/${genome_version}.masked \
    --genomeFastaFiles ${genome_version}.mask/${genome_version}.fa.masked \
    --runThreadN 21

#mkdir /scratch/fbeaudry/alignments
for sample in 243_liver_ACTTGA_L003 243_liver_ACTTGA_L005 251_heart_GATCAG_L003 251_heart_GATCAG_L005 251_kidney_TTAGGC_L003 251_kidney_TTAGGC_L005 253_heart_GCCAAT_L003 253_heart_GCCAAT_L005 253_kidney_ACAGTG_L003 253_kidney_ACAGTG_L005 253_ovary_CGATGT_L003 253_ovary_CGATGT_L005 
do
    STAR --genomeDir ${genome_version}.masked \
        --readFilesIn RNAseq/NC_FSJ_${sample}_R1_val_1.fq.gz RNAseq/NC_FSJ_${sample}_R2_val_2.fq.gz \
        --twopassMode Basic --outFileNamePrefix alignments/${sample}. \
        --runThreadN 21 --readFilesCommand zcat

    java -jar picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT \
        I=alignments/${sample}.Aligned.out.sam \
        O=alignments/${sample}.tmp.add.bam  \
        RGID=${sample} RGLB=${sample} RGPL=RNA RGPU=unit1 RGSM=${sample}

    samtools sort -n alignments/${sample}.tmp.add.bam >alignments/${sample}.sort.bam

done

module unload perl/5.24.0/b1
module load braker/2.1.6 diamond/0.9.24/b2 prothint/2.6.0 genemark/b1 augustus/3.4.0 perl/5.32.0/b1

mkdir BRAKER1/${run}

perl braker.pl \
    --gff3 --cores=21 --etpmode \
    --genome=${genome_version}.fa.masked \
    --DIAMOND_PATH=diamond/0.9.24/b2/bin/ \
    --CDBTOOLS_PATH=cdbfasta/ \
    --AUGUSTUS_CONFIG_PATH=augustus/Augustus/config/ \
    --PROTHINT_PATH=prothint/2.6.0/bin/ \
    --GENEMARK_PATH=genemark/4.68 \
    --AUGUSTUS_SCRIPTS_PATH=augustus/3.4.0/scripts \
    --AUGUSTUS_BIN_PATH=augustus/3.4.0/bin \
    --workingdir=BRAKER1/${run} \
    --bam=alignments/243_liver_ACTTGA_L003.sort.bam,alignments/253_heart_GCCAAT_L003.sort.bam,alignments/243_liver_ACTTGA_L005.sort.bam,alignments/253_heart_GCCAAT_L005.sort.bam,alignments/251_heart_GATCAG_L003.sort.bam,alignments/253_kidney_ACAGTG_L003.sort.bam,alignments/251_heart_GATCAG_L005.sort.bam,alignments/253_kidney_ACAGTG_L005.sort.bam,alignments/251_kidney_TTAGGC_L003.sort.bam,alignments/253_ovary_CGATGT_L003.sort.bam,alignments/251_kidney_TTAGGC_L005.sort.bam,alignments/253_ovary_CGATGT_L005.sort.bam \

module unload python3/3.8.3
module load busco/5.2.2 bedtools

cp BRAKER1/${run}/braker.gtf genome_files/${genome_version}.BRAKER1.gtf

awk '$3 ~ "gene" {print $1"\t"$2"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' BRAKER1/${run}/braker.gtf > BRAKER1/${run}/braker.${run}.gtf

bedtools getfasta -fi ${genome_version}.mask/${genome_version}.fa.masked \
    -bed BRAKER1/${run}/braker.${run}.gtf \
    -nameOnly -fo BRAKER1/${run}/braker.${run}.gtf.fa

busco -m transcriptome -i BRAKER1/${run}/braker.${run}.gtf.fa \
    -o ${genome_version}.BRAKER1_busco.tsv  -l aves_odb10 -f -c 21

##STRINGTIE

#Stringtie steps by Mitch mitch.lok3@gmail.com 

#SBATCH -c 21 --mem=50G
module unload python3/3.8.3
module load busco/5.2.2 bedtools

awk '$3 ~ "transcript" {print $1"\t"$2"\t"$12"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ${genome_version}.stringtieMerged.gtf |  awk 'split($3,a,"\"") {print $1"\t"$2"\t"a[2]"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' >${genome_version}.stringtieMerged.adj.gtf

bedtools getfasta -fi ${genome_version}.mask/${genome_version}.fa.masked \
    -bed ${genome_version}.stringtieMerged.adj.gtf -nameOnly \
    -fo ${genome_version}.stringtie.fa

busco -m transcriptome -i ${genome_version}.stringtie.fa -o ${genome_version}.stringtie.busco -l aves_odb10 -f -c 21



