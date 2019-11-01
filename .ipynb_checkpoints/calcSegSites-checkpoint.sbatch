#!/bin/bash
#
#SBATCH -J calcSegSites # A single job name for the array
#SBATCH -p shared # Partition
#SBATCH -n 1 # cores
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # Running time of 0 days 02 hrs 00 min
#SBATCH --mem 8000 # Memory request of 8 GB
#SBATCH -o tmp/segSites_%A_%a.out # Standard output (what is written to console)
#SBATCH -e tmp/segSites_%A_%a.err # Standard error (errors written to console)
#SBATCH --account=debivort_lab
#SBATCH --open-mode=append

#load some modules
module load Anaconda3/5.0.1-fasrc02

#grab your python environment
source activate py_env

#get your population
NAME_FILES=($(ls -1 sample_names/herit_pops))
OUT_NAME=$(echo ${NAME_FILES[$SLURM_ARRAY_TASK_ID]} | egrep -o '.._herit|.._var')

#run the python script
python ~/Seq-Data/bootstrap_indvs.py sample_names/herit_pops/ca_herit_files_v2.txt sample_names/herit_pops/${NAME_FILES[$SLURM_ARRAY_TASK_ID]} $OUT_NAME