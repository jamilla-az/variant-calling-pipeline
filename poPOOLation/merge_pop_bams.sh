#!/bin/bash
#
#SBATCH -J bam_pop_merge # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-0:30 # Running time of 30 min
#SBATCH --mem 6000 # Memory request of 6 GB
#SBATCH -o /n/debivort_lab/Jamilla_seq/tmp/bamMerge_%j.out # Standard output (what is written to console)
#SBATCH -e /n/debivort_lab/Jamilla_seq/tmp/bamMerge_%j.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load module
module load samtools
module load python/3.6.3-fasrc02

DIR_PATH=/n/debivort_lab/Jamilla_seq/poPOOLation/bam_lists

#find unique lines and pick n randomly - do this in python to save headache
N=5
python ~/Seq-Data/poPOOLation/sample_unique_lines.py $DIR_PATH/va_bam_herit.txt $N $DIR_PATH/va_bam_herit_merge.txt

#merge BAMs from indviduals in those lines into a 'population'
samtools merge -b va_bam_herit_merge.txt 