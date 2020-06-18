#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 09:22:08 2020

this script produces sliding window replicates within chromosome-boundries of phylip-files.
that can be used for HyDe runs

@author: Thomsn
"""
import os
import pandas as pd

CURDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/hyde'
WINDOW_SIZE = 30

os.chdir(CURDIR)

snp_file = '../2020_04_vitis_snps_summary.csv'
data_file = '2020_04_vitis_USA-data.txt'

snp_sum = pd.read_csv(snp_file)
snp_sum = snp_sum.sort_values(by=['chromosome','locus'])
repdir = f'window_replicates_size{WINDOW_SIZE}'
if not os.path.exists(repdir):
    os.mkdir(repdir)
    os.mkdir(f'{repdir}/data')
    
total_frame = pd.DataFrame(columns = ['header', 'sequence'])
with open(data_file, 'r') as daf:
    soup = daf.readlines()
    for i, molecule in enumerate(soup):
        separ = molecule.split(' ')
        head = ' '.join(separ[:-1])
        sequence = separ[-1].split('\n')[0]
        total_frame.loc[str(i),:] = [head, sequence]
        
# for chro in set(snp_sum['chromosome']):
#     sel_loci = snp_sum.copy()[snp_sum['chromosome'] == chro]
#     for coord in range(len(sel_loci)+1-WINDOW_SIZE):
#         sel_window = sel_loci.iloc[coord:coord+WINDOW_SIZE,:].copy()
#         wdw_names = list(sel_window['locus'].values)
#         out_name = f'{repdir}/data/wdsize{WINDOW_SIZE}_chr{chro}_{wdw_names[0]}_{wdw_names[-1]}-data.txt'
#         if not os.path.exists(repdir):
#             with open(out_name, 'w') as ousi:
#                 for _, row in total_frame.iterrows():
#                     head = row['header']
#                     sequence = row['sequence']
#                     split_seq = ''.join([base for i, base in enumerate(sequence) if i in sel_window['Unnamed: 0'].values])
#                     ousi.write(f'{head} {split_seq}\n')
        
snp_names = [f'{snn["chromosome"]}_{snn["locus"]}' for _, snn in snp_sum.iterrows()]
snp_set = pd.DataFrame(columns=snp_names)
for _, row in total_frame.iterrows():
    head = row['header']
    sequence = row['sequence']
    sequencesplit = list(sequence)
    split_seq = [sequencesplit[i] for i in snp_sum['Unnamed: 0'].values]
    snp_set.loc[head,:] = split_seq
    # split_seq = ''.join([base for i, base in enumerate(sequence) if i in snp_sum['Unnamed: 0'].values])
    

for chro in set(snp_sum['chromosome']):
    sel_loci = snp_sum.copy()[snp_sum['chromosome'] == chro]
    for coord in range(len(sel_loci)+1-WINDOW_SIZE):
        sel_window = sel_loci.iloc[coord:coord+WINDOW_SIZE,:].copy()
        wdw_names = list(sel_window['locus'].values)
        out_name = f'{repdir}/data/wdsize{WINDOW_SIZE}_chr{chro}_{wdw_names[0]}_{wdw_names[-1]}-data.txt'
        selsnp_set = snp_set.copy().loc[:,list(snp_sum['chromosome'] == chro)]
        print(selsnp_set)
        if not os.path.exists(out_name):
            with open(out_name, 'w') as ousi:
                for _, row in selsnp_set.iterrows():
                    head = row.name
                    split_seq = ''.join(row[coord:coord+WINDOW_SIZE])
                    # split_seq = ''.join([base for i, base in enumerate(sequence) if i in sel_window['Unnamed: 0'].values])
                    ousi.write(f'{head} {split_seq}\n')
        




