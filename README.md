# Variant Calling Pipeline

#### Scripts for BWA-MEM, picard, GATK, vcftools to use with whole-genome variant calling

Explicit details of the whole process can be found in NovaSeq Analysis Pipeline.md

**alignment** contains the scripts for aligning the FASTQ files to the reference genome and processing the BAM files using BWA-MEM and picard tools. 

**genotyping** contains the scripts for calling genotypes and generating the final VCF file using GATK tools. 

**variant_processing** contains custom Python, vcftools, and snpEff scripts to process the VCF files. 

**sample_scripts** contains test scripts for working with single samples only (used for debugging tools and estimating memory/time usage). 

