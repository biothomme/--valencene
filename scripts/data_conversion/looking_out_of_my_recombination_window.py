#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 08:52:35 2020

this script should give an overview about non-recombinant regions within
the vitis dataset

@author: Thomsn
"""
import os
import pandas as pd

def get_all_snps(wdw_frame):
    return [snp for subl in wdw_frame['snp_names'] for snp in subl]

def find_mutual_rwdws(frame1, frame2):
    next_snps = get_all_snps(frame2)
    current_snps = get_all_snps(frame1)
    mutual_snps = [snp for snp in next_snps if snp in current_snps]
    rec_frames = frame1['snp_names']
    csel_wdws = [co for co, rec in enumerate(rec_frames) if all([True for r in rec if r in mutual_snps])]
    old_wdws = frame1.iloc[csel_wdws,]
    rec_frames = frame2['snp_names']
    dsel_wdws = [co for co, rec in enumerate(rec_frames) if all([True for r in rec if r in mutual_snps])]
    new_wdws = frame2.iloc[dsel_wdws,]
    new_half = []
    old_half = []
    for wdw in old_wdws['snp_names']:
        new_pdt = [i for i, npd in enumerate(new_wdws['snp_names']) if sum([True for el in npd if el in wdw]) == len(npd)]
        for npdt in new_pdt:
            # print(f'{wdw} --- {new_wdws["snp_names"].values[npdt]}')
            new_half.append(npdt)
    for wdw in new_wdws['snp_names']:
        old_pdt = [i for i, opd in enumerate(old_wdws['snp_names']) if (sum([True for el in wdw if el in opd]) == len(opd) & len(opd) < len(wdw))]
        for opdt in old_pdt:
                old_half.append(opdt)
    nww = new_wdws.iloc[new_half,:].copy()
    oww = old_wdws.iloc[old_half,:].copy()
    sum([1 for nnnn in nww['snp_names'] if nnnn in list(oww['snp_names'].values)])
    new_recs = pd.concat([nww,oww])
    new_recs['chrb_start'] = pd.to_numeric(new_recs['chrb_start'].copy())
    new_recs['chromosome'] = [f'0{rest}' if len(rest) == 1 else rest for rest in new_recs['chromosome'].copy()]
    out_recs = new_recs.sort_values(by=['chromosome', 'chrb_start']).copy()
    return out_recs


BIG_DIR = 'Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/non_recomb_data/'
os.chdir(BIG_DIR)
HEADER = ['chromosome',
          'chrb_start',
          'chrb_end',
          'length',
          'snps_count',
          'snp_names']

tax_list = [diro for diro in os.listdir() if '.det' in diro]
tax_list = tax_list[:2] + tax_list[-2:]
total_frame = pd.Series()
for tax_file in tax_list:
    with open(tax_file, 'r') as tf:
        soup = tf.readlines()
        window_frame = pd.DataFrame(columns = HEADER)
        for coco, molecule in enumerate(soup):
            molecule = molecule.replace('\n', '')
            parts = [atom for atom in molecule.split(' ') if len(atom) > 0]
            if coco > 0:
                window_frame.loc[str(coco),:] = parts
        window_frame['snp_names'] = [nm.split('|') for nm in window_frame['snp_names']]
    total_frame[tax_file] = window_frame

rec_fr = total_frame[0].copy()
for count in range(1,total_frame.size):
    print(count)
    rec_fr = find_mutual_rwdws(rec_fr, total_frame[count]).copy()
rec_fr.to_csv('max_rec_wdws_vvinifera.csv')
