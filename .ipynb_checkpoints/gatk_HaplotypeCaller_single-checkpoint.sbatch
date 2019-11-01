#!/bin/bash
#
#SBATCH -J haplo_call # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 20 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 0-4:00 # Running time of 4 hrs 00 min
#SBATCH --mem 16000 # Memory request of 16 GB
#SBATCH -o tmp/haplo_%j.out # Standard output (what is written to console)
#SBATCH -e tmp/haplo_%j.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

cd 03_testGvcfs

#load some modules
module load jdk/1.8.0_45-fasrc01 #java
module load gatk/4.0.2.1-fasrc01 #GATK

#make some variables to store paths
BAM_DIR=~/Seq-Data/02_testBams
REF_FILE=~/Seq-Data/00_genome/dmel-all-chromosome-r6.28.fasta

#run haplotype caller - can specify chr with -L option
gatk HaplotypeCaller \
    --java-options "-Xmx16g -XX:ParallelGCThreads=20" \
    -R $REF_FILE \
	-I $BAM_DIR/BERK_1_1217_S218_final.bam \
	-O BERK_1_1217_S218.2L.20.raw.g.vcf \
    -L 2L \
    --emit-ref-confidence GVCF \
    --min-pruning 1 \
    --min-dangling-branch-length 1 
