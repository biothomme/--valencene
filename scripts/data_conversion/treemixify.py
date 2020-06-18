#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:25:15 2020

this script should convert a biallelic matrix and population 
assignments to get an output that can be fed into treemix

@author: Thomsn
"""

import os
import pandas as pd

BIGDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
SNP_FILE = 'sorted_snp_set.txt'
POP_FILE = '2020_04_vitis_USA_2_summary.csv'

os.chdir(BIGDIR)

snp_array = pd.read_csv(SNP_FILE)
snp_array.index = list(snp_array.iloc[:,0].values)
snp_array = snp_array.drop(columns='Unnamed: 0')

pop_array = pd.read_csv(POP_FILE)
pop_array.index = list(pop_array.iloc[:,0].values)
pop_array = pop_array.drop(columns='Unnamed: 0')
cutnms = lambda x: x.lower().replace(' ','').replace('.','')
pop_array['Pop-celine'] = [cutnms(x) for x in pop_array['Pop-celine']]

with open('vitis_treemix.txt', 'w') as vt:
    header = []
    for pop in set(pop_array['Pop-celine'].values):
        header.append(pop)
    vt.write(f'{" ".join(header)}\n')
        
    for _, line in snp_array.iterrows():
        pop_counts = []
        for pop in set(pop_array['Pop-celine'].values):
            inds = list(pop_array[pop_array['Pop-celine'].values == pop]['cut_names'].values)
            total_inds = len(inds) * 2
            alleles = line[inds]
            nonal = sum(alleles)
            domal = total_inds - nonal
            pop_counts.append(f'{domal},{nonal}')
        vt.write(f'{" ".join(pop_counts)}\n')
        
import numpy as np
[int(i**(2**.5))+1 for i in range(20)]
