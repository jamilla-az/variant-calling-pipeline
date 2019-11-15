#!/usr/bin/python -tt

#import os
import sys
import argparse
import random
import pandas as pd
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description='Bootstrap estimator for number of segregating sites')

    parser.add_argument('base_file', type=str,
                        help='table of base (CA for herit, BERK for var) individuals labels')
    parser.add_argument('indv_file', type=str,
                        help='table of X pop individuals labels')
    parser.add_argument('output_label', type=str,
                        help='pop label for output file')
    parser.add_argument(
        'var_type', type=str, help='variant type label for output file: all, genic, exon, nonsyn')

    args = parser.parse_args()
    print(f'using the following list of individuals: {args.indv_file}')

    # dataframe of individuals as baseline for bootstrapping
    base_df = pd.read_csv(args.base_file, sep='\t', header=None,
                          index_col=0, names=['loc', 'line', 'num'])

    # dataframe of individuals to use from X pop
    indv_df = pd.read_csv(args.indv_file, sep='\t', header=None,
                          index_col=0, names=['loc', 'line', 'num'])

    # dataframe of genotypes x individuals
    geno_df = pd.read_csv(
        f'/n/debivort_lab/Jamilla_seq/final_vcfs/output_files/genotype_matrix_{args.var_type}.table', sep='\t')
    # print(geno_df.head())

    # find the individuals with called genotypes in the list of individuals to use
    # this are the individuals we will bootstrap from
    boot_indv = [i for i in indv_df.index if i in geno_df.columns]
    indv_df = indv_df.loc[boot_indv]
    ca_indv = [i for i in base_df.index if i in geno_df.columns]
    seg_sites = []  # set empty list to store bootstrap output

    # 100 bootstraps
    random.seed(0)
    for n in np.arange(0, 100):
        print(f'sampling individuals with replacement: iter {n+1}')
        # bootstrap your base individuals and count up the number of individuals in each line
        ids = random.choices(ca_indv, k=len(ca_indv))
        lines = base_df['line'].loc[ids].value_counts().to_frame(
            name='num_indv')

        # bootstrap your other individuals to follow the same line/num_indvs breakdown
        inds_list = []
        l = random.sample(list(indv_df['line'].unique()), len(lines))
        for k in np.arange(0, len(lines)):
            inds = random.choices(
                indv_df.loc[indv_df['line'] == l[k]].index, k=lines['num_indv'].values[k])
            inds_list = inds_list + inds

        # filter your genotype x indv matrix to include only bootstrapped indvs
        temp_df = geno_df[inds_list]
        ss = 0  # set seg sites counter to 0

        print('tallying genotypes')
        geno_counts = temp_df.apply(pd.Series.value_counts, axis=1).fillna(0)
        geno_counts.columns = ['hom-ref', 'het', 'hom-alt']
        geno_counts = geno_counts.astype('int')

        print('calculating number of seg sites')
        # calculating number of segregating sites
        for i in np.arange(0, len(geno_counts)):
            # num ref alleles = num hom ref genotypes
            hom_ref = geno_counts['hom-ref'].iloc[i]
            # num alt alleles = num hom alt genotypes
            hom_alt = geno_counts['hom-alt'].iloc[i]
            # random ref/alt sampling for hets
            het_sample = np.random.choice(2, geno_counts['het'].iloc[i])
            hom_ref += np.sum(het_sample == 0)  # add het ref alleles
            hom_alt += np.sum(het_sample == 1)  # add het alt alleles

            if (hom_ref > 0) & (hom_alt > 0):  # if that add a segregating site
                ss += 1

        print('appending to list')
        seg_sites.append(ss)

    print(
        f'saving list as output_files/boot_{args.output_label}_segsites_{args.var_type}.table')
    #print(pd.DataFrame({'seg_sites': seg_sites}))
    pd.DataFrame({'seg_sites': seg_sites}).to_csv(
        f'output_files/boot_{args.output_label}_segsites_{args.var_type}.table', sep='\t', index=False)


if __name__ == '__main__':
    main()
