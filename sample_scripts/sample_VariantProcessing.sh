#!/bin/bash

cd 05_testVCF

#load some modules
module load jdk/1.8.0_45-fasrc01 #java

#how many SNPs have ReadPosRankSum as a field
grep 'ReadPosRankSum' wildFlies_chr2_test.vcf | grep -v '#' | wc -l

#half of SNPs do not have this field, so can't filter by it well
#this is because the two samples are both homozygous at this SNP! 
grep -v 'ReadPosRankSum' wildFlies_chr2_test.vcf | grep -v '#' | wc -l


#how many do not have SNPs have QD as a field (34)
grep -v 'QD' wildFlies_chr2_test.vcf | grep -v '#' | wc -l

#show me a sample of those lines
grep -v 'QD' wildFlies_chr2_test.vcf | grep -v '##' > no_QD.vcf

#hard filter variants based on quality metrics
~/gatk-4.1.3.0/gatk VariantFiltration \
-R ../00_genome/dmel-all-chromosome-r6.28.fasta \
-V wildFlies_chr2_test.vcf \
-O wildFlies_chr2_test_filtered.vcf \
--filter-expression "QD < 2.0 || MQ < 40.0 || SOR > 3.0 || FS > 60.0" \
--filter-name "bad"


#how many SNPs passed the filter
grep -o 'PASS' wildFlies_chr2_test_filtered.vcf | wc -l

#how many SNPs failed the filter
grep -o 'bad' wildFlies_chr2_test_filtered.vcf | wc -l

#how many do not have SNPs have QD as a field - still 34, and are marked with . in the Filter field 
grep -v 'QD' wildFlies_chr2_test_filtered.vcf | grep -v '#' | wc -l

#sample of passing SNPs 
grep 'PASS' wildFlies_chr2_test_filtered.vcf | head -n 10 > sample_filtered.vcf

## let's check out vcftools

module load vcftools
module load vcflib
module load bcftools

#calculate mean depth per individual and return an idepth file
vcftools --vcf wildFlies_chr2_test_filtered.vcf --depth --out meanDepth
vcftools --vcf wildFlies_chr2_test_filtered.vcf

#get VCF statistics
bcftools stats wildFlies_chr2_test_filtered.vcf > vcfStats.txt

#missing data per individual
vcftools --vcf wildFlies_chr2_test_filtered.vcf --remove-filtered-all --minDP 2 --out missingData --missing-indv

#relatedness
vcftools --vcf wildFlies_chr2_test_filtered.vcf --remove-filtered-all --minDP 2 --out relatedScore --relatedness2q