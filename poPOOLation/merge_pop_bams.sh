#!/bin/bash
#
#SBATCH -J bam_pop_merge # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-1:00 # Running time of 1 hrs 00 min
#SBATCH --mem 6000 # Memory request of 6 GB
#SBATCH -o /n/debivort_lab/Jamilla_seq/tmp/bamMerge_%j.out # Standard output (what is written to console)
#SBATCH -e /n/debivort_lab/Jamilla_seq/tmp/bamMerge_%j.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

# $1 = name of pop e.g. ca, ma, fl
# $2 = expt subset e.g. herit or var
# $3 = num lines to subsample 

#load module
module load samtools
module load python/3.6.3-fasrc02
module load Anaconda3/5.0.1-fasrc02
module load jdk/1.8.0_172-fasrc01

DIR_PATH=/n/debivort_lab/Jamilla_seq/poPOOLation/bam_lists
OUT_PATH=/n/debivort_lab/Jamilla_seq/poPOOLation/pop_bams

#find unique lines and pick n randomly - do this in python to save headache
N=$3
python ~/Seq-Data/poPOOLation/sample_unique_lines.py $DIR_PATH/$1_bam_$2.txt $N $DIR_PATH/$1_bam_$2_merge.txt

#merge BAMs from indviduals in those lines into a 'population'
samtools merge -b $DIR_PATH/$1_bam_$2_merge.txt $OUT_PATH/$1_$2_merged.bam

#sort it and output index file and sorted BAM
java -jar ~/picard.jar SortSam \
I=$OUT_PATH/$1_$2_merged.bam \
O=$OUT_PATH/$1_$2_merged_sorted.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true