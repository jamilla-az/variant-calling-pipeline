#!/bin/bash

#HaplotypeCaller will make your BAM files into GVCFs (an intermediate and slow step) - this helps if you are doing this in portions. You can make a whole bunch of GVCFs and then analyze them all together vs. having to re-do the whole thing everytime you add new data if you go straight to a VCF

#we will do this on file BERK_1_1217_S218_final.bam in 02_testBams

#get into your directory
cd 03_testGvcfs

#load some modules
module load jdk/1.8.0_45-fasrc01 #java
module load gatk/4.0.2.1-fasrc01 #GATK

#make some variables to store paths
BAM_DIR=~/Seq-Data/02_testBams
REF_FILE=~/Seq-Data/00_genome/dmel-all-chromosome-r6.28.fasta

#run haplotype caller - can specify chr with -L option
gatk HaplotypeCaller \
    --java-options "-Xmx8g -XX:ParallelGCThreads=2" \
    -R $REF_FILE \
	-I $BAM_DIR/BERK_1_1217_S218_final.bam \
	-O BERK_1_1217_S218.2L.raw.g.vcf \
    -L 2L \
    --emit-ref-confidence GVCF \
    --min-pruning 1 \
    --min-dangling-branch-length 1 
