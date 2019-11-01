#!/bin/bash
#
#SBATCH -J bwa_alignRerun # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-2:00 # Running time of 1 hrs 00 min
#SBATCH --mem 8000 # Memory request of 8 GB
#SBATCH -o tmp/bwaRerun_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/bwaRerun_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append


#sample file names: BERK_1_1217_S218_L001_R1_001.fastq.gz, BERK_1_1217_S218_L001_R2_001.fastq.gz
#flowcell barcode: HFFLJDRXX
#sample code: bwa mem -M -t 1 -R '@RG\tID:{FLOWCELL}.{LANE}\tPU:{FLOWCELL_BARCODE}.{LANE}.{SAMPLE}\tSM:{SAMPLE}\tPL:{PLATFORM}\tLB{LIBRARY}' <genome_prefix> <reads_1.fq> <reads_2.fq> > <samplename_bwa.sam>

cd 02_bams

#get list of fastq files - unique sample names only (FROM FAILED SAMPLES ONLY)
FILES=($(ls -1 /n/scratchlfs/debivort_lab/Jamilla/01_fastqs/failed_samples | sed 's/_R.*//' | uniq))

#grab file of interest
FILE_NAME=${FILES[$SLURM_ARRAY_TASK_ID]}

#extract everything before the lane number as sample
SAMPLE=$(echo $FILE_NAME | sed 's/_L00.*//')

#extract lane number
LANE=$(echo $FILE_NAME | grep -o 'L00.' | tail -c 2)

### fix fastqs to remove empty reads (reads < 15bp) ###
module load jdk

java -jar ~/Trimmomatic-0.35/trimmomatic-0.35.jar PE ../01_fastqs/failed_samples/${SAMPLE}_L00${LANE}_R1_001.fastq.gz ../01_fastqs/failed_samples/${SAMPLE}_L00${LANE}_R2_001.fastq.gz ../01_fastqs/failed_samples/${SAMPLE}_L00${LANE}_R1_001_paired.fastq.gz ../01_fastqs/failed_samples/${SAMPLE}_L00${LANE}_R1_001_unpaired.fastq.gz ../01_fastqs/failed_samples/${SAMPLE}_L00${LANE}_R2_001_paired.fastq.gz ../01_fastqs/failed_samples/${SAMPLE}_L00${LANE}_R2_001_unpaired.fastq.gz MINLEN:15


#align using bwa mem 

module load bwa

bwa mem -M -t 1 -R "@RG\tID:HFFLJDRXX.$LANE\tPU:HFFLJDRXX.$LANE.$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA\tLB:$SAMPLE" \
../00_genome/dmel_r6 \
../01_fastqs/failed_samples/${SAMPLE}_L00${LANE}_R1_001_paired.fastq.gz \
../01_fastqs/failed_samples/${SAMPLE}_L00${LANE}_R2_001_paired.fastq.gz > ${SAMPLE}_${LANE}.sam

# convert to bam

#sort by coordinate and create an index 
java -jar ~/picard.jar SortSam \
I=${SAMPLE}_${LANE}.sam \
O=${SAMPLE}_${LANE}_sorted.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true

#collect alignment metrics

#java -jar ~/picard.jar CollectAlignmentSummaryMetrics \
I=${SAMPLE}_${LANE}_sorted.bam \
R=../00_genome/dmel-all-chromosome-r6.28.fasta \
METRIC_ACCUMULATION_LEVEL=SAMPLE \
METRIC_ACCUMULATION_LEVEL=READ_GROUP \
O=${SAMPLE}_${LANE}.alignment_metrics.txt


