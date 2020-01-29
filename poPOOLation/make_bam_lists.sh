#grab list of BAM files that come from a specific population - heritability

ls -1 /n/debivort_lab/Jamilla_seq/final_bams/*.bam | grep 'CA_[0-9]*_F' > /n/debivort_lab/Jamilla_seq/poPOOLation/bam_lists/ca_bam_herit.txt

#grab list of BAM files that come from a specific population - variability

ls -1 /n/debivort_lab/Jamilla_seq/final_bams/*.bam | grep -e 'FL_[0-9]*_[0-9]\|MIA_[0-9]*_[0-9]' > /n/debivort_lab/Jamilla_seq/poPOOLation/bam_lists/fl_bam_var.txt