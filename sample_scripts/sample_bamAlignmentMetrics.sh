#!/bin/bash
#sample processing of 1 bam file - found in 02_testBams

cd 02_testBams

#step 1: sort and index your .sam file
#step 2: collect summary alignment metrics

#load in java module

module load jdk

#specify your file and pull out sample name only
FILES=($(ls -1 *.bam))
SAMPLE=$( echo ${FILES[1]} | sed 's/_sorted.bam//')

#collect alignment metrics

java -jar ~/picard.jar CollectAlignmentSummaryMetrics \
I=${FILES[1]} \
R=../00_genome/dmel-all-chromosome-r6.28.fasta \
METRIC_ACCUMULATION_LEVEL=SAMPLE \
METRIC_ACCUMULATION_LEVEL=READ_GROUP \
O=${SAMPLE}.alignment_metrics.txt