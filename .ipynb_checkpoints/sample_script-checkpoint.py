#!/usr/bin/python -tt

#import os
#import sys
import pandas as pd
import numpy as np
#import random

def main():
    #load table
    var_df = pd.read_csv('output_files/chr2L_flHerit_fullGeno.table',sep='\t').drop(['Unnamed: 9'], axis =1)
    #print(len(var_df))
    allele_array = np.empty([len(var_df), 2], dtype=int) #initialize empty array
    seg_sites = []
    #100 bootstraps
    for n in np.arange(0,100):
        ss = 0
        #iterate over each site
        for i in np.arange(0, len(var_df)):
            allele_array[i,0] = var_df['HOM-REF'].iloc[i] # num ref alleles = num hom ref genotypes
            allele_array[i,1] = var_df['HOM-VAR'].iloc[i] # num alt alleles = num hom alt genotypes
            het_sample = np.random.choice(2, var_df['HET'].iloc[i]) #random ref/alt sampling for hets
            allele_array[i,0] += np.sum(het_sample == 0) #add het ref alleles
            allele_array[i,1] += np.sum(het_sample == 1) #add het alt alleles

            if (allele_array[i,0] > 0) & (allele_array[i,1] > 0):
                ss += 1
                
        seg_sites.append(ss)
    
    var_df = var_df.assign(REF = allele_array[:,0],
                 ALT = allele_array[:,1]) #assign new cols
    
    var_df.to_csv('output_files/chr2L_flHerit_fullGeno_updated.table', sep ='\t', index = False)
    
    pd.DataFrame({'num_seg_sites':seg_sites}).to_csv('output_files/chr2L_flherit_segsites.table', sep='\t', index = False)
    

if __name__ == '__main__':
    main()
    