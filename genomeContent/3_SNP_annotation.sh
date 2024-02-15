##annotate SNPs

module load plink/1.9 java/1.8.0_111

# build annotation, using file output from geneCombo.R
java -Xmx8g -jar snpEff.jar build -gtf22 -v ${genome_version}

# convert SNP file from plink to vcf
plink --noweb --allow-no-sex --dog --nonfounders --file geno.anon --recode vcf --out geno.anon

# format input file
Rscript --vanilla convertVCF.R

# run snpEff
java -Xmx8g -jar snpEff.jar FSJ.v2_1 geno.anon.faux.vcf > geno.anon.faux.ann.vcf

# reformat 
awk '$1 !~ "#" {print $3"\t"$8"\t"$1"\t"$2}' geno.anon.faux.ann.vcf | awk 'split($2,a,"|") {print $1"\t"a[2]"\t"a[3]"\t"a[5]"\t"$3"\t"$4}' > geno.anon.snp_annotation.txt

# in R, combine with information in snp.pos.txt and save as snp_annotation.txt