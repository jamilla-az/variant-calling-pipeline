#!/bin/bash
#sample processing of 1 bam file - found in 02_testBams

cd 02_testBams

#step 1: sort and index your .sam file

#load in java module

module load jdk

#specify your file and pull out sample name only
FILES=($(ls -1 *.sam))
SAMPLE=$( echo ${FILES[1]} | sed 's/\.sam//')

#sort by coordinate and create an index 
java -jar ~/picard.jar SortSam \
I=${SAMPLE}.sam \
O=${SAMPLE}_sorted.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true
