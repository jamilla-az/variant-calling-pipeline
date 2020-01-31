#!/bin/bash
#
#SBATCH -J bam_sort # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-0:30 # Running time of 30 min
#SBATCH --mem 6000 # Memory request of 6 GB
#SBATCH -o /n/debivort_lab/Jamilla_seq/tmp/bamSort_%A_%a.out # Standard output (what is written to console)
#SBATCH -e /n/debivort_lab/Jamilla_seq/tmp/bamSort_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load java module
module load Java/1.8

#pick subsampled BAM file
FILES=($(ls -1 /n/debivort_lab/Jamilla_seq/subsampled_bams/*.bam))
SAMPLE=${FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(echo $SAMPLE | sed 's/_subsample.bam//' | sed 's/\/n\/debivort_lab\/Jamilla_seq\/subsampled_bams\///')

#sort it and output index file and sorted BAM
java -jar ~/picard.jar SortSam \
I=${SAMPLE} \
O=/n/debivort_labs/Jamilla_seq/subsampled_bams/sorted_bams/${SAMPLE_NAME}_sampled_sort.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true