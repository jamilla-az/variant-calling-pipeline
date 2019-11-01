#!/bin/bash

#debug any failed sam to bam conversions post-alignment

#go into the sam/bam directory in scratchlfs
cd /n/scratchlfs/debivort_lab/Jamilla/02_bams

#check that all fastq files were made into sam files
#/n/scratchlfs/debivort_lab/Jamilla/01_fastqs

#get list of fastq files - unique sample_lane names only
ls -1 /n/scratchlfs/debivort_lab/Jamilla/01_fastqs/*.fastq.gz | sed 's/_R.*//' | uniq > fastqs.txt
#get list of sam files - unique sample_lane names only
ls -1 *.sam | sed 's/.\.sam/L00&/' | sed 's/_*.sam//'> sams.txt
#get list of bam files - unique sample_lane names only
ls -1 *_sorted.bam | sed 's/._sorted\.bam/L00&/' | sed 's/.sorted\.bam//'> bams.txt

#find fastqs that do not match a sam file (alignment failed)
grep -vf sams.txt fastqs.txt > fastqs_nosam.txt

#find sam files that do not match a bam file (sam to bam conversion failed, no output)
grep -vf bams.txt sams.txt > sam_nobam.txt

#pull out fastqs without a sam file for later realignment
NOSAM_FILES=($(ls -1 ../01_fastqs | grep -f fastqs_nosam.txt))

#pull out sam files without a bam file
sed 's/L00//' sam_nobam.txt > test.txt
ls -1 *.sam | sed 's/\.sam//' > allsam_files.txt
NOBAM_FILES=($(grep -f allsam_files.txt test.txt))

#load java module
module load jdk

#validate each sam file without a bam file
mkdir bam_debug

for FILE in ${NOBAM_FILES[@]}
do
    #SAMPLE=$(echo $FILE | sed 's/_*.sam//')
    java -jar ~/picard.jar ValidateSamFile \
    I=$FILE.bam \
    O=bam_debug/$FILE.validate.txt \
    MODE=VERBOSE
done

#find empty bam files
find *.bam -type f -empty > empty_bams.txt

#pull out sam files with an empty bam file
sed 's/_sorted.*//' empty_bams.txt > test.txt
ls -1 *.sam | sed 's/\.sam//' > allsam_files.txt
EMPTYBAM_FILES=($(grep -f allsam_files.txt test.txt))


### search for validation files with errors ###

wc -l 02_bams/bam_debug/* | grep -o '0'  #search for empty validation files

grep -vRl 'No errors found' ../02_bams/bam_debug/



#pull out fastq files with an empty bam file
sed 's/._sorted\.bam/L00&/' empty_bams.txt | sed 's/.sorted\.bam//'> test.txt
sed 's/\/n\/scratchlfs\/debivort_lab\/Jamilla\/01_fastqs\///' fastqs.txt | uniq > test_fastqs.txt
EMPTY_FASTQS=($(grep -f test.txt test_fastqs.txt))

#copy them to the failed samples directory
for FILE in ${EMPTY_FASTQS[@]}
do
    cp -n ../01_fastqs/${FILE}* ../01_fastqs/failed_samples
done


### fix fastqs to remove empty reads ###
java -jar ~/Trimmomatic-0.35/trimmomatic-0.35.jar PE ~/failed_samples/MA_8_F3_1045_S46_L002_R1_001.fastq.gz ~/failed_samples/MA_8_F3_1045_S46_L001_R1_001.fastq.gz ~/failed_samples/MA_8_F3_1045_S46_L001_R2_001_paired.fastq.gz ~/failed_samples/MA_8_F3_1045_S46_L001_R2_001_unpaired.fastq.gz ~/failed_samples/MA_8_F3_1045_S46_L001_R1_001_paired.fastq.gz ~/failed_samples/MA_8_F3_1045_S46_L001_R1_001_unpaired.fastq.gz MINLEN:15

### realign ###
module load bwa

bwa mem -M -t 1 -R "@RG\tID:HFFLJDRXX.1\tPU:HFFLJDRXX.1.MA_8_F3_1045_S46\tSM:MA_8_F3_1045_S46\tPL:ILLUMINA\tLB:MA_8_F3_1045_S46" \
~/failed_samples/00_genome/dmel_r6 \
~/failed_samples/MA_8_F3_1045_S46_L001_R1_001_paired.fastq.gz \
~/failed_samples/MA_8_F3_1045_S46_L001_R2_001_paired.fastq.gz > test.sam


#sort by coordinate and create an index 
module load jdk 

java -jar ~/picard.jar SortSam \
I=test.sam \
O=test_sorted.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true

### testing #####
java -jar ~/picard.jar ValidateSamFile \
I=test_sorted.bam \
O=test.validate.txt \
MODE=SUMMARY

### search for validation files with errors ###

#wc -l 02_bams/final_bams/valid_files/* | grep -o '0'  #search for empty validation files
grep -vRl 'No errors found' ../02_bams/final_bams/valid_files/

