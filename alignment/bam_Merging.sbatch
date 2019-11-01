#!/bin/bash
#
#SBATCH -J bam_merge # A single job name for the array
#SBATCH -p serial_requeue # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-0:30 # Running time of 30 min
#SBATCH --mem 6000 # Memory request of 6 GB
#SBATCH -o tmp/bamMerge_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/bamMerge_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

cd 02_bams

#load in java module

module load jdk

#specify your file and pull out sample name only
FILES=($(ls -1 *.bam | sed 's/_[12]_sorted.bam//' | uniq))
SAMPLE=${FILES[$SLURM_ARRAY_TASK_ID]}

#merge bam files from the sample and output to separate directory called merged_bams in 02_bams

java -jar ~/picard.jar MergeSamFiles \
I=${SAMPLE}_1_sorted.bam \
I=${SAMPLE}_2_sorted.bam \
O=merged_bams/${SAMPLE}_merged.bam \
CREATE_INDEX=true

