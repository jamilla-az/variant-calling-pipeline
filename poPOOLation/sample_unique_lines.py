#!/usr/bin/python -tt

import sys
import argparse
import pandas as pd
import numpy as np
import random


def main():
    parser = argparse.ArgumentParser(
        description='Sample N random lines & output list of individual BAMs')

    parser.add_argument('bam_list', type=str,
                        help='list of BAM files')
    parser.add_argument('num_lines', type=int,
                        help='number of lines to sample')
    parser.add_argument('output_file', type=str,
                        help='file path + name for output txt file')

    args = parser.parse_args()
    print(f'using the following list of individuals: {args.bam_list}')

    # find the unique lines
    lines = pd.read_csv(args.bam_list, sep='\/',
                        header=None, engine='python')
    lines = lines.iloc[:, 6].str.split("_", n=2, expand=True)
    uniq_lines = lines.iloc[:, 1].unique()

    # sample n unique lines
    print(f'sampling this many lines: {args.num_lines}')
    sampled_lines = random.sample(list(uniq_lines), args.num_lines)

    # find indices of BAM files that belong to these lines
    file_idx = []
    for line in sampled_lines:
        file_idx.extend(lines[1].loc[lines[1] == line].index)

    bam_files = pd.read_csv(args.bam_list, sep='\t',
                            header=None, engine='python')
    selected_bams = bam_files.loc[file_idx]

    # write selected BAM files to a txt file
    selected_bams.to_csv(args.output_file, index=False, header=False, sep='\t')


if __name__ == '__main__':
    main()
