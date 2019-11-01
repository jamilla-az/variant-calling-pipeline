#!/usr/bin/python -tt

from cyvcf2 import VCF
import pandas as pd
import numpy as np


def main():
    print('starting code')
    vcf = VCF('allChr_fullGeno.recode.vcf', gts012=True)
    
    geno_mat = np.empty([46988, len(vcf.samples)], dtype=int) #initialize empty array
    chr_list = []
    pos_list = []
    i = 0
    for variant in vcf:
        chr_list.append(variant.CHROM)
        pos_list.append(variant.start)
        geno_mat[i,:] = variant.gt_types #HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UKNOWN=3
        i +=1 
    
    geno_df = pd.DataFrame(geno_mat, columns = vcf.samples)
    geno_df.assign(chr = chr_list, pos = pos_list).to_csv('output_files/genotype_matrix.table', sep='\t', index = False)
    
    print(pos_list[0])
    print('new file made')
    
if __name__ == '__main__':
    main()


