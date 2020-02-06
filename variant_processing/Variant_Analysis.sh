#!/bin/bash

#load some modules
module load jdk/1.8.0_45-fasrc01 #java
module load vcftools
#module load Anaconda3/5.0.1-fasrc02 #for python packages

#gather VCFs for whole genome
java -jar ~/picard.jar GatherVcfs \
I=wildFlies_2L.vcf.gz \
I=wildFlies_2R.vcf.gz \
I=wildFlies_3L.vcf.gz \
I=wildFlies_3R.vcf.gz \
I=wildFlies_4.vcf.gz \
I=wildFlies_X.vcf.gz \
O=wildFlies_all.vcf.gz

#index file
~/gatk-4.1.3.0/gatk IndexFeatureFile \
-F wildFlies_all.vcf.gz

#hard filter variants based on quality metrics
~/gatk-4.1.3.0/gatk VariantFiltration \
-R ../00_genome/dmel-all-chromosome-r6.28.fasta \
-V wildFlies_all.vcf.gz \
-O wildFlies_all_filtered.vcf.gz \
--filter-expression "QD < 2.0 || MQ < 40.0 || SOR > 3.0 || FS > 60.0" \
--filter-name "bad"

#calc missingness per indv
vcftools --gzvcf wildFlies_all_filtered.vcf.gz \
--remove-filtered-all \
--missing-indv \
--out output_files/missing_indv_all

#calculate mean depth per individual and return an idepth file
vcftools --gzvcf wildFlies_all_filtered.vcf.gz \
--remove-filtered-all \
--depth \
--out output_files/meanDepth_all_filt

#filter out individuals with very low coverage 
awk '{if ($3 > 2) print $1}' output_files/meanDepth_all_filt.idepth | grep -v 'INDV' > sample_names/filt_samples_all.txt

#filter out all sites that have any missing data and keep original info field
vcftools --gzvcf wildFlies_all_filtered.vcf.gz \
--remove-filtered-all \
--keep sample_names/filt_samples_all.txt \
--recode \
--recode-INFO-all \
--max-missing 1 \
--out allChr_fullGeno

#filter out all sites that have any missing data and keep original info field - biallelic sites only
vcftools --gzvcf wildFlies_all_filtered.vcf.gz \
--remove-filtered-all \
--keep sample_names/filt_samples_all.txt \
--recode \
--recode-INFO-all \
--max-missing 0.5 \
--remove-indels \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.05 \
--out allChr_halfGeno_biallelic

#filter for a particular population and chr
vcftools --gzvcf wildFlies_3R.vcf.gz \
--keep sample_names/herit_pops/va_herit_files_v2.txt \
--recode \
--recode-INFO-all \
--out va_herit_3R

#calc missingness per indv for VA flies on chr 3R
vcftools --vcf va_herit_3R.recode.vcf \
--missing-indv \
--out output_files/missing_indv_va_X_4_nonfilt

#calc mean depth per indv for VA flies on chr 3R
vcftools --vcf va_herit_3R.recode.vcf \
--depth \
--out output_files/meanDepth_va_3R_nonfilt




#make sample name files

#MA flies - heritability only
grep '^MA_[0-9]*_F[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > ma_herit_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/herit_pops/ma_herit_files.txt > sample_names/herit_pops/ma_herit_files_v2.txt

#MA flies - variability only
grep '^MA_[0-9]*_[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > ma_var_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/var_pops/ma_var_files.txt > sample_names/var_pops/ma_var_files_v2.txt

#FL flies - heritability only
grep '^FL_[0-9]*_F[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > fl_herit_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/herit_pops/fl_herit_files.txt > sample_names/herit_pops/fl_herit_files_v2.txt

#FL flies - variability only
grep -E '^(FL|MIA)_[0-9]*_[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > fl_var_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/var_pops/fl_var_files.txt > sample_names/var_pops/fl_var_files_v2.txt

#VA flies - heritability only
grep '^VA_[0-9]*_F[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > va_herit_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/herit_pops/va_herit_files.txt > sample_names/herit_pops/va_herit_files_v2.txt

#VA flies - variability only
grep '^CM[0-9]*' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > va_var_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/var_pops/va_var_files.txt > sample_names/var_pops/va_var_files_v2.txt

#TX flies - heritability only
grep '^TX_[0-9]*_F[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > tx_herit_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/herit_pops/tx_herit_files.txt > sample_names/herit_pops/tx_herit_files_v2.txt

#TX flies - variability only
grep '^TX_[0-9]*_[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > tx_var_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/var_pops/tx_var_files.txt > sample_names/var_pops/tx_var_files_v2.txt

#PA flies - heritability only
grep '^PA_[0-9]*_F[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > pa_herit_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/herit_pops/pa_herit_files.txt > sample_names/herit_pops/pa_herit_files_v2.txt

#PA flies - variability only
grep '^PA_[0-9]*_[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > pa_var_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/var_pops/pa_var_files.txt > sample_names/var_pops/pa_var_files_v2.txt

#CA flies - heritability only
grep '^CA_[0-9]*_F[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > ca_herit_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/herit_pops/ca_herit_files.txt > sample_names/herit_pops/ca_herit_files_v2.txt

#CA flies - variability only
grep '^CA_[0-9]*_[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > ca_var_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/var_pops/ca_var_files.txt > sample_names/var_pops/ca_var_files_v2.txt

#BERK flies - variability only
grep '^BERK_[0-9]*_[0-9]' ../gvcf_files.txt | sed 's/.raw.g.vcf//' | uniq > berk_var_files.txt
awk 'BEGIN {FS ="_";OFS ="\t"};{print $0,$1,$2,$3}' sample_names/var_pops/berk_var_files.txt > sample_names/var_pops/berk_var_files_v2.txt
