#!/bin/bash
#
#SBATCH -J haplo_call # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 20 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 0-18:00 # Running time of 18 hrs 00 min
#SBATCH --mem 16000 # Memory request of 16 GB
#SBATCH -o tmp/haplo_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/haplo_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load some modules
module load jdk/1.8.0_45-fasrc01 #java

#make some variables to store paths
BAM_DIR=02_bams/final_bams
REF_FILE=00_genome/dmel-all-chromosome-r6.28.fasta

#get list of final BAM files
FILES=($(ls -1 $BAM_DIR/*.bam | sed 's/02_bams\/final_bams\///'))

#grab file of interest
FILE_NAME=${FILES[$SLURM_ARRAY_TASK_ID]}

#extract sample name
SAMPLE=$(echo $FILE_NAME | sed 's/_final.bam//')


#run haplotype caller - whole genome, can specify chr with -L option
~/gatk-4.1.3.0/gatk HaplotypeCaller \
    --java-options "-Xmx16g -XX:ParallelGCThreads=20" \
    -R $REF_FILE \
	-I $BAM_DIR/$FILE_NAME \
	-O 03_gcvfs/${SAMPLE}.raw.g.vcf \
    --emit-ref-confidence GVCF \
    --min-pruning 1 \
    --min-dangling-branch-length 1 
