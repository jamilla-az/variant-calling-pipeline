#!/bin/bash
#
#SBATCH -J pileup # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-4:00 # Running time of 4hrs 00 min
#SBATCH --mem 6000 # Memory request of 6 GB
#SBATCH -o /n/debivort_lab/Jamilla_seq/tmp/pileup_%A_%a.out # Standard output (what is written to console)
#SBATCH -e /n/debivort_lab/Jamilla_seq/tmp/pileup_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

IN_PATH=~/popoolation-code
OUT_PATH=/n/debivort_lab/Jamilla_seq/poPOOLation/pop_bams

#pick your pop
LOC_ARRAY=(va ma fl ca tx pa)
POP=herit
LOC=${LOC_ARRAY[$SLURM_ARRAY_TASK_ID]}

#load module
module load samtools
module load perl/5.26.1-fasrc01

#convert population bam into pileup file
samtools mpileup -f ~/Seq-Data/00_genome/dmel-all-chromosome-r6.28.fasta -q 20 -o $OUT_PATH/${LOC}_${POP}.pileup $OUT_PATH/${LOC}_${POP}_merged_sorted.bam 

#mask + filter out SNPs near indels
perl $IN_PATH/basic-pipeline/identify-genomic-indel-regions.pl --input $OUT_PATH/${LOC}_${POP}.pileup --output $OUT_PATH/${LOC}_${POP}_indels.gtf --indel-window 5
perl $IN_PATH/basic-pipeline/filter-pileup-by-gtf.pl --gtf $OUT_PATH/${LOC}_${POP}_indels.gtf --input $OUT_PATH/${LOC}_${POP}.pileup --output $OUT_PATH/${LOC}_${POP}_filt.pileup

#calculate theta with a sliding window of 50kb
POOL_SIZE=$(wc -l ../bam_lists/${LOC}_bam_${POP}_merge.txt | grep -o '[0-9]*')

perl $IN_PATH/Variance-sliding.pl \
--input ${LOC}_${POP}_filt.pileup --output ${LOC}_${POP}.theta \
--measure theta --window-size 50000 --step-size 50000 \
--min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size $POOL_SIZE