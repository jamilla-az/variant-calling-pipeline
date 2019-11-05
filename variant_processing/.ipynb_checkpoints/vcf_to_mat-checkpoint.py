#!/usr/bin/python -tt

from cyvcf2 import VCF
import argparse
import pandas as pd
import numpy as np


def main():
    print('parsing arguments')
    parser = argparse.ArgumentParser(description='Convert VCF to matrix of genotypes')
    
    parser.add_argument('num_vcf_lines', type=int, help="num lines in VCF file (grep -v '#' <VCF> | wc -l)")
    parser.add_argument('VCF_file', type=str, help='VCF to convert to matrix')
    parser.add_argument('out_suffix', type=str, help='file suffix to add to output file: all, genic, exon, nonsyn')
    
    args = parser.parse_args()
    print(f'using VCF {args.VCF_file} with {args.num_vcf_lines}')                   
    print('starting code')
    vcf = VCF(args.VCF_file, gts012=True)
    
    geno_mat = np.empty([args.num_vcf_lines, len(vcf.samples)], dtype=int) #initialize empty array
    chr_list = []
    pos_list = []
    i = 0
    for variant in vcf:
        chr_list.append(variant.CHROM)
        pos_list.append(variant.start)
        geno_mat[i,:] = variant.gt_types #HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UKNOWN=3
        i +=1 
    
    geno_df = pd.DataFrame(geno_mat, columns = vcf.samples)
    geno_df.assign(chr = chr_list, pos = pos_list).to_csv(f'output_files/genotype_matrix_{args.out_suffix}.table', sep='\t', index = False)
    
    #print(pos_list[0])
    print('new file made')
    
if __name__ == '__main__':
    main()


