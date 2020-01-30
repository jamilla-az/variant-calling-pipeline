#!/bin/bash
#
#SBATCH -J bam_downsample # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-0:15 # Running time of 15 min
#SBATCH --mem 6000 # Memory request of 6 GB
#SBATCH -o tmp/bamSample_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/bamSample_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load in samtools module
module load samtools

#specify the minimum read count
MIN_READS=100000

#specify bam file
FILES=($(ls -1 /n/debivort_lab/Jamilla_seq/final_bams/*.bam))
SAMPLE=${FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(echo $SAMPLE | sed 's/_final.bam//')

#count the num of reads in the current bam
CURR_READS=$(samtools view -c $SAMPLE)

#get fraction to downsample to 
n=$(bc -l <<< "scale=2;$MIN_READS/$CURR_READS")

#downsample and save new bam
samtools view -bs 42$n $SAMPLE > /n/debivort_lab/Jamilla_seq/subsampled_bams/${SAMPLE_NAME}_subsample.bam
