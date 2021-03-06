#!/bin/bash
#
#SBATCH -J bam_downsample # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-0:30 # Running time of 30 min
#SBATCH --mem 6000 # Memory request of 6 GB
#SBATCH -o /n/debivort_lab/Jamilla_seq/tmp/bamSample_%A_%a.out # Standard output (what is written to console)
#SBATCH -e /n/debivort_lab/Jamilla_seq/tmp/bamSample_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load in samtools module
module load samtools

#specify the minimum read count
MIN_READS=$(samtools view -c /n/debivort_lab/Jamilla_seq/final_bams/CA_3_F1_1004_S5_final.bam)

#specify bam file
#FILES=($(ls -1 /n/debivort_lab/Jamilla_seq/final_bams/*.bam))
#SAMPLE=${FILES[$SLURM_ARRAY_TASK_ID]}

#specify bam files to rerun
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/Seq-Data/poPOOLation/bams_to_sample.txt)
SAMPLE_NAME=$(echo $SAMPLE | sed 's/_final.bam//' | sed 's/\/n\/debivort_lab\/Jamilla_seq\/final_bams\///')

#count the num of reads in the current bam
CURR_READS=$(samtools view -c $SAMPLE)

#get fraction to downsample to 
n=$(bc -l <<< "scale=2;$MIN_READS/$CURR_READS")

#downsample and save new bam
samtools view -bs 42$n $SAMPLE > /n/debivort_lab/Jamilla_seq/subsampled_bams/${SAMPLE_NAME}_subsample.bam
