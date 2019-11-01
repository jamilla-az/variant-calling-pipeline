#!/bin/bash

#will make the mapping between sample and the GVCF file name for GenomicsDBImport. 
#From the GATK tool description: "The sample map is a tab-delimited text file with sample_name--tab--path_to_sample_vcf per line. Using a sample map saves the tool from having to download the GVCF headers in order to determine the sample names. Sample names in the sample name map file may have non-tab whitespace, but may not begin or end with whitespace."

#collect all GVCF file names into a text file
cd 03_gvcfs #in scratchlfs Jamilla directory

ls -1 *.raw.g.vcf > ../gvcf_files.txt #store in Jamilla directory

#pull out sample names using awk and attach relative file path to GVCF file 
cd ..
awk -F. '{print $1"\t../03_gvcfs/"$0}' gvcf_files.txt > test.txt 