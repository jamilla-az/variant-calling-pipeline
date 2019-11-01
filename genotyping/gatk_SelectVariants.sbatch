#!/bin/bash
#
#SBATCH -J selVars # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 0-00:30 # Running time of 0 days 00 hrs 30 min
#SBATCH --mem 12000 # Memory request of 12 GB
#SBATCH -o tmp/selVars_%j.out # Standard output (what is written to console)
#SBATCH -e tmp/selVars_%j.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load some modules
module load jdk/1.8.0_45-fasrc01 #java

#pick your chromosome
#CHR=(2L 2R 3L 3R X 4)
#SAMPLE=${CHR[$SLURM_ARRAY_TASK_ID]}

cd 05_testVCF

~/gatk-4.1.3.0/gatk SelectVariants \
--java-options "-Xmx8g" \
-R ../00_genome/dmel-all-chromosome-r6.28.fasta \
-V gendb://../04_genDB/my_database_2R \
-L 2R:1323000-1333000 \
-O wildFlies_2R.vcf.gz
