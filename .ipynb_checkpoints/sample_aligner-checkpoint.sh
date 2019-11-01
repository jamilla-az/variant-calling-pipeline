#!/bin/bash
#sample data alignment script using data in 01_testFastqs directory - 1 sample

#sample file names: BERK_1_1217_S218_L001_R1_001.fastq.gz, BERK_1_1217_S218_L001_R2_001.fastq.gz
#flowcell barcode: HFFLJDRXX
#sample code: bwa mem -M -t 1 -R '@RG\tID:{FLOWCELL}.{LANE}\tPU:{FLOWCELL_BARCODE}.{LANE}.{SAMPLE}\tSM:{SAMPLE}\tPL:{PLATFORM}\tLB{LIBRARY}' <genome_prefix> <reads_1.fq> <reads_2.fq> > <samplename_bwa.sam>

#set the directory
cd ~/Seq-Data

#grab file name of interest
FILE_NAME='BERK_1_1217_S218_L001_R1_001.fastq.gz'

#extract everything before the lane number as sample
SAMPLE=$(echo $FILE_NAME | sed 's/_L00.*//')

#extract lane number
LANE=$(echo $FILE_NAME | grep -o 'L00.' | tail -c 2)

#align using bwa mem 
cd 02_testBams

module load bwa

bwa mem -M -t 1 -R "@RG\tID:HFFLJDRXX.$LANE\tPU:HFFLJDRXX.$LANE.$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA\tLB:$SAMPLE" \
../00_genome/dmel_r6 \
../01_testFastqs/${SAMPLE}_L00${LANE}_R1_001.fastq.gz \
../01_testFastqs/${SAMPLE}_L00${LANE}_R2_001.fastq.gz > ${SAMPLE}_${LANE}.sam


