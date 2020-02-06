#!/bin/bash
#
#SBATCH -J genDB # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 0-8:00 # Running time of 8 hrs 00 min
#SBATCH --mem 12000 # Memory request of 12 GB
#SBATCH -o tmp/genDB_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/genDB_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load some modules
module load jdk/1.8.0_45-fasrc01 #java

#pick your chromosome
CHR=(2L 2R 3L 3R X 4)
SAMPLE=${CHR[$SLURM_ARRAY_TASK_ID]}

cd 04_genDB #make sure there is a tmp dir there

~/gatk-4.1.3.0/gatk GenomicsDBImport \
--java-options "-Xmx8g" \
--genomicsdb-workspace-path my_database_${SAMPLE} \
--batch-size 50 \
-L $SAMPLE \
--tmp-dir=$SCRATCH/debivort_lab/tmp \
--sample-name-map ../gvcf.sample_map
