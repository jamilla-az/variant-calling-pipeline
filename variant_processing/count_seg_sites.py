#!/usr/bin/python -tt

# import os
import sys
import argparse
import pandas as pd
import numpy as np
import random


def main():
    parser = argparse.ArgumentParser(
        description='Estimator for number of segregating sites')

    parser.add_argument('indv_file', type=str,
                        help='table of X pop individuals labels')
    parser.add_argument('output_label', type=str,
                        help='pop label for output file')
    parser.add_argument(
        'var_type', type=str, help='variant type label for output file: all, genic, exon, nonsyn')

    args = parser.parse_args()
    print(f'using the following list of individuals: {args.indv_file}')

    # dataframe of individuals to use from X pop
    indv_df = pd.read_csv(args.indv_file, sep='\t', header=None,
                          index_col=0, names=['loc', 'line', 'num'])

    # dataframe of genotypes x individuals
    geno_df = pd.read_csv(
        f'/n/debivort_lab/Jamilla_seq/final_vcfs/output_files/genotype_matrix_{args.var_type}.table', sep='\t')

    # pull out the unique lines
    lines = indv_df['line'].value_counts()
    # remove any lines that don't have 4 indiv labels recorded
    lines_to_use = lines.loc[lines == 4].index
    seg_sites = []  # initialize empty list for storing num of seg sites
    lines_used = []  # initialize empty list for storing lines used

    random.seed(0)
    for line in lines_to_use:
        print(f"pulling out individuals from line {line}")
        # find the individuals with called genotypes in the list of individuals to use
        inds_list = [i for i in indv_df.loc[indv_df["line"]
                                            == line].index if i in geno_df.columns]
        print(inds_list)
        if len(inds_list) == 4:
            # use line if has 4 individuals with called geno
            lines_used.append(line)
            # filter your genotype x indv matrix to include only those indvs
            temp_df = geno_df[inds_list]
            ss = 0  # set seg sites counter to 0

            print('tallying genotypes')
            geno_counts = temp_df.apply(
                pd.Series.value_counts, axis=1).fillna(0)
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

            print('appending num seg sites to list')
            seg_sites.append(ss)

        else:
            print(
                f'Skipping...line {line} has < 4 ({len(inds_list)}) individuals')
    print(
        f'saving list as output_files/count_{args.output_label}_segsites_{args.var_type}.table')

    pd.DataFrame({'loc': [args.output_label]*len(seg_sites), 'lines': lines_used, 'seg_sites': seg_sites}).to_csv(
        f'output_files/count_{args.output_label}_segsites_{args.var_type}.table', sep='\t', index=False)


if __name__ == '__main__':
    main()
