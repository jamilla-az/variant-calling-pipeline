#!/bin/bash
#sample processing of 1 bam file - found in 02_testBams

cd 02_testBams

#step 1: sort and index your .sam file
#step 2: collect summary alignment metrics
#step 3: merge BAMs from the sample sample

#load in java module

module load jdk

#specify your file and pull out sample name only
FILES=($(ls -1 *.bam))
SAMPLE=$( echo ${FILES[0]} | sed 's/_[12]_sorted.bam//')

#merge bam files from the sample sample
java -jar ~/picard.jar MergeSamFiles \
I=${SAMPLE}_1_sorted.bam \
I=${SAMPLE}_2_sorted.bam \
O=${SAMPLE}_merged.bam \
CREATE_INDEX=true