#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 08:44:14 2020

this script is made to summarize the output of each of the subsequent floating
window replicates. Furthermore it can be used to plot ABBA against BABA patterns.

@author: Thomsn
"""

__author__ = 'thomas m. huber'
__email__ = ['thomas.huber@evobio.eu']

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# obtain list of abbas and babas of my dataset
def mybabaisabba(rep_list):
    rep_names = [rp.split(f'wdsize{SNP_WDW_SIZE}_')[1].split('-data')[0] for rp in rep_list]
    abbas = pd.DataFrame(columns = rep_names)
    babas = pd.DataFrame(columns = rep_names)
    aabbs = pd.DataFrame(columns = rep_names)
    for wdw_name, wdw_file in zip(rep_names, rep_list):
        abba_list = []
        baba_list = []
        aabb_list = []
        with open(wdw_file, 'r') as wf:
            soup = wf.readlines()
            for i, molecule in enumerate(soup):
                if i != 0:
                    atoms = molecule.split('\t')
                    total_val = [float(at.split('e+')[0]) * 10**float(at.split('e+')[0]) if 'e+' in at else float(int(at.split('.')[0])) for at in atoms[6:]]
                    abbab = total_val[8]
                    babab = total_val[6]
                    # babab = ( total_val[6] + total_val[5] )/2
                    aabbb = total_val[3]
                    total_sum = sum(total_val)
                    abba = abbab / total_sum
                    baba = babab / total_sum
                    aabb = aabbb / total_sum
                    abba_list.append(abba)
                    baba_list.append(baba)
                    aabb_list.append(aabb)
        abbas[wdw_name] = abba_list
        babas[wdw_name] = baba_list
        aabbs[wdw_name] = aabb_list
    return abbas, babas, aabbs

SNP_WDW_SIZE = 30
BIGDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/hyde/admixed_asia/wdw_reps_size30/out_vbg267'
os.chdir(BIGDIR)

rep_list = [diro for diro in os.listdir() if '-ind-bonf.txt' in diro]
# rep_list.sort()
rep_names = [rp.split(f'wdsize{SNP_WDW_SIZE}_')[1].split('-data')[0] for rp in rep_list]

startbase = [int(nm.split('_')[1]) for nm in rep_names]
endbase = [int(nm.split('_')[2]) for nm in rep_names]
chroms = [int(nm.split('_')[0].split('chr')[1]) for nm in rep_names]
rep_frame = pd.DataFrame({'file': rep_list, 
                          'start': startbase, 
                          'end': endbase, 
                          'chr': chroms})
rep_frame = rep_frame.sort_values(by = ['chr','start'])


# abbas, babas = mybabaisabba(rep_frame['file'].values)
abbas, babas, aabbs = mybabaisabba(rep_frame['file'].values)


for chrom in set(list(rep_frame['chr'].values)):
    selector = list(rep_frame['chr'] == chrom)
    print(abbas.loc[:,selector].mean().values[0:10])
    y = rep_frame.copy()[selector]['start'].values
    fig = plt.figure(figsize=(30, 8))
    ax = plt.subplot(111)
    # for (_, abba), (_, baba) in zip(abbas.loc[:,selector].iterrows(), babas.loc[:,selector].iterrows()):
    for (_, abba), (_, baba), (_, aabb) in zip(abbas.loc[:,selector].iterrows(), 
                                               babas.loc[:,selector].iterrows(),
                                               aabbs.loc[:,selector].iterrows()):
        ax.plot(y, abba.values, color = 'teal', alpha = .15)
        ax.plot(y, baba.values, color = 'sienna', alpha = .15)
        ax.plot(y, aabb.values, color = 'limegreen', alpha = .15)
    # ax.plot(y, abbas.loc[:,selector].mean().values, color = 'teal', label = 'V. v. east (ABBA)', linewidth = 3)
    # ax.plot(y, babas.loc[:,selector].mean().values, color = 'sienna', label = 'V. sylv. west (BABA)', linewidth = 3)
    # ax.plot(y, babas.loc[:,selector].mean().values, color = 'sienna', label = 'V. sylv. west (BABA) + V sylv. (BABB)', linewidth = 3)
    ax.plot(y, abbas.loc[:,selector].mean().values, color = 'teal', label = 'V. piasezkii (BABA)', linewidth = 3)
    ax.plot(y, babas.loc[:,selector].mean().values, color = 'sienna', label = 'V. acerifolia (ABBA)', linewidth = 3)
    ax.plot(y, aabbs.loc[:,selector].mean().values, color = 'limegreen', label = 'V. sylv. west (AABB)', linewidth = 3)
    ax.legend(loc='upper left', title='frequency of shared private sites with ...')
#    ax.set_title(f'signs of admixture into V. v. west along chromosome {chrom} (Outg.: V. sylv east).')
#    ax.set_title(f'signs of admixture into Vbg12 along chromosome {chrom}.')
    ax.set_xlabel(f'base count (of first SNP of {SNP_WDW_SIZE} SNP window)')
    ax.set_ylabel('frequency of site pattern from total site patterns')
    plt.savefig(f'{chrom}_plot.pdf')

 #   list(abbas.loc[:,selector].mean().values)
