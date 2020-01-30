#! /bin/bash

#use to move labeled low coverage indvs to a different folder 
num_indvs=$(wc -l low_coverage_indvs.txt | grep -o '[0-9]*')
a=$((num_indvs+1)) #add 1 for proper indexing

for n in $(seq 1 $a)
do
    file=$(sed "${n}q;d" low_coverage_indvs.txt)
    mv /n/debivort_lab/Jamilla_seq/final_bams/${file}_final.bai /n/debivort_lab/Jamilla_seq/low_coverage_bams
done