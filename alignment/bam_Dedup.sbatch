#!/bin/bash
#
#SBATCH -J bam_dedup # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-0:30 # Running time of 30 min
#SBATCH --mem 6000 # Memory request of 6 GB
#SBATCH -o tmp/bamDedup_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/bamDedup_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

cd 02_bams/merged_bams

#load in java module

module load jdk

#specify your file and pull out sample name only

FILES=($(ls -1 *_merged.bam))
FILE_NAME=${FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$( echo $FILE_NAME | sed 's/_merged.bam//')

#mark duplicates

java -jar ~/picard.jar MarkDuplicates \
TMP_DIR=tmp \
I=$FILE_NAME \
O=${SAMPLE}_dedup.bam \
METRICS_FILE=${SAMPLE}.dedup.metrics.txt \
REMOVE_DUPLICATES=false \
TAGGING_POLICY=All

#sort again and output to final_bams directory in 02_bams

java -jar ~/picard.jar SortSam \
I=${SAMPLE}_dedup.bam \
O=../final_bams/${SAMPLE}_final.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true

#validate your BAM file and output to final_bams directory in 02_bams

java -jar ~/picard.jar ValidateSamFile \
I=../final_bams/${SAMPLE}_final.bam \
O=../final_bams/${SAMPLE}.validate.txt \
MODE=SUMMARY