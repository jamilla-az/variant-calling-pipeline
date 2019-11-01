#!/bin/bash

#load some modules
module load jdk/1.8.0_45-fasrc01 #java

#pick your chromosome
CHR='chr2'

cd 05_testVCF

~/gatk-4.1.3.0/gatk GenotypeGVCFs \
--java-options "-Xmx6g" \
-R ../00_genome/dmel-all-chromosome-r6.28.fasta \
-V gendb://../04_testGenDB/my_database \
-O wildFlies_${CHR}_test.vcf.gz \
--tmp-dir=tmp