#!/bin/bash
#sample processing of 1 bam file - found in 02_testBams

cd 02_testBams

#step 1: sort and index your .sam file
#step 2: collect summary alignment metrics
#step 3: merge BAMs from the sample sample
#step 4: mark duplicates and re-sort

#load in java module

module load jdk

#specify your file and pull out sample name only
FILES=($(ls -1 *_merged.bam))
SAMPLE=$( echo ${FILES[0]} | sed 's/_merged.bam//')

#mark duplicates
java -jar ~/picard.jar MarkDuplicates \
TMP_DIR=tmp \
I=${FILES[0]} \
O=${SAMPLE}_dedup.bam \
METRICS_FILE=${SAMPLE}.dedup.metrics.txt \
REMOVE_DUPLICATES=false \
TAGGING_POLICY=All

#sort again
java -jar ~/picard.jar SortSam \
I=${SAMPLE}_dedup.bam \
O=${SAMPLE}_final.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true