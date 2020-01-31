#grab list of BAM files that come from a specific population - heritability

ls -1 /n/debivort_lab/Jamilla_seq/subsampled_bams/*.bam | grep 'TX_[0-9]*_F' > /n/debivort_lab/Jamilla_seq/poPOOLation/bam_lists/tx_bam_herit.txt

#grab list of BAM files that come from a specific population - variability

ls -1 /n/debivort_lab/Jamilla_seq/subsampled_bams/*.bam | grep -e 'BERK_[0-9]*_[0-9]' > /n/debivort_lab/Jamilla_seq/poPOOLation/bam_lists/berk_bam_var.txt