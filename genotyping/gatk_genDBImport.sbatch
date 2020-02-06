#!/bin/bash
#
#SBATCH -J genDB # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 2-00:00 # Running time of 2 days 00 hrs 00 min
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

#genDB database must be created in POSIX compliant filesystem
#stat -f -c %T . to check current directory - lustre is compliant, but nfs is not

OUT_PATH=/n/debivort_lab/Jamilla_seq

IN_PATH=$SCRATCH/debivort_lab #make sure there is a tmp dir here and in the genDB folder 

~/gatk-4.1.3.0/gatk GenomicsDBImport \
--java-options "-Xmx8g" \
--genomicsdb-workspace-path $IN_PATH/genDB/my_database_${SAMPLE} \
--batch-size 50 \
-L $SAMPLE \
--tmp-dir=$SCRATCH/debivort_lab/tmp \
--sample-name-map $OUT_PATH/gvcf.sample_map
