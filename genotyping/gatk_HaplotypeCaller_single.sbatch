#!/bin/bash
#
#SBATCH -J haplo_call # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 2 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 0-20:00 # Running time of 20 hrs 00 min
#SBATCH --mem 16000 # Memory request of 16 GB
#SBATCH -o tmp/haplo_%j.out # Standard output (what is written to console)
#SBATCH -e tmp/haplo_%j.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

cd 03_testGvcfs

#load some modules
module load jdk/1.8.0_45-fasrc01 #java
#module load gatk/4.0.2.1-fasrc01 #GATK

#make some variables to store paths
BAM_DIR=~/Seq-Data/02_testBams
REF_FILE=~/Seq-Data/00_genome/dmel-all-chromosome-r6.28.fasta

#run haplotype caller - can specify chr with -L option
~/gatk-4.1.3.0/gatk HaplotypeCaller \
--java-options "-Xmx16g -XX:ParallelGCThreads=2" \
-R $REF_FILE \
-I $BAM_DIR/BERK_1_1217_S218_final.bam \
-O BERK_1_1217_S218.all.2.raw.g.vcf \
--emit-ref-confidence GVCF \
--min-pruning 1 \
--min-dangling-branch-length 1 
