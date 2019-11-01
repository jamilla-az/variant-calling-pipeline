#!/bin/bash
#
#SBATCH -J genoVCF # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 4-00:00 # Running time of 4 days 00 hrs 00 min
#SBATCH --mem 12000 # Memory request of 12 GB
#SBATCH -o tmp/genoVCF_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/genoVCF_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load some modules
module load jdk/1.8.0_45-fasrc01 #java

#pick your chromosome
#chr2R,3L,3R corrupted genDB
CHR=(2L 2R 3L 3R X 4)
SAMPLE=${CHR[$SLURM_ARRAY_TASK_ID]}

cd 05_vcf

~/gatk-4.1.3.0/gatk GenotypeGVCFs \
--java-options "-Xmx8g" \
-R ../00_genome/dmel-all-chromosome-r6.28.fasta \
-V gendb://../04_genDB/my_database_${SAMPLE} \
-O wildFlies_${SAMPLE}.vcf.gz \
--tmp-dir=tmp
