#!/bin/bash
#sample processing of 1 bam file - found in 02_testBams

cd 02_testBams

#step 1: sort and index your .sam file
#step 2: collect summary alignment metrics
#step 3: merge BAMs from the sample sample
#step 4: mark duplicates and re-sort
#step 5: validate your BAM file

#load in java module

module load jdk

#specify your file and pull out sample name only
FILES=($(ls -1 *_final.bam))
SAMPLE=$( echo ${FILES[0]} | sed 's/_final.bam//')

java -jar ~/picard.jar ValidateSamFile \
I=${FILES[0]} \
O=${SAMPLE}.validate.txt \
MODE=SUMMARY