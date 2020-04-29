#!/bin/bash
#
#SBATCH -J calcSegSites # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # Running time of 0 days 01 hrs 00 min
#SBATCH --mem 8000 # Memory request of 8 GB
#SBATCH -o tmp/segSites_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/segSites_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load some modules
module load Anaconda3/5.0.1-fasrc02

#grab your python environment
source activate py_env

DIR_PATH=/n/debivort_lab/Jamilla_seq/final_vcfs

#get your population
NAME_FILES=($(ls -1 $DIR_PATH/sample_names/herit_pops))
OUT_NAME=$(echo ${NAME_FILES[$SLURM_ARRAY_TASK_ID]} | egrep -o '.._herit|.._var|...._var')
VAR_TYPE='all'

#run the python script
# var - var_pops/berk_var_files_v2.txt
# herit - herit_pops/ca_herit_files_v2.txt

python ~/Seq-Data/variant_processing/bootstrap_indvs.py $DIR_PATH/sample_names/herit_pops/ca_herit_files_v2.txt $DIR_PATH/sample_names/herit_pops/${NAME_FILES[$SLURM_ARRAY_TASK_ID]} $OUT_NAME $VAR_TYPE
