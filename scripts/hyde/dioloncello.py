#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 13:14:30 2020

@author: Thomsn
"""

import os
import pandas as pd
from itertools import combinations
import argparse
import numpy as np


# os.chdir('/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/hyde')
# takes infile and mapfile and coverts it to pandas dataframes
def to_pededeff(infile,
                head=True):
    with open(infile) as table:
        table_lines = table.readlines()
        for i, line in enumerate(table_lines):
            spline = line.split('\t')
            spline = [field.split('\n')[0] for field in spline]
            if i == 0:
                if head:
                    header = spline
                    d_dat = pd.DataFrame(index = header)
                else:
                    header = ['ind', 'pop']
                    d_dat = pd.DataFrame(data=spline, index = header)
            else:
                frame = pd.DataFrame(data=spline, index=header)
                d_dat = pd.concat([d_dat,frame], axis=1)
    return d_dat

# calculates ABBA-BABA statistics
def mama_mia(all_dat): 
    dstat = (all_dat.loc['ABBA'] - all_dat.loc['ABAB']) /\
    (all_dat.loc['ABBA'] + all_dat.loc['ABAB'])
    
    return dstat

# calculates Hils statistic (Kubatko and Chifman 2019)
def up_the_hils(all_dat): 
    total_inv = sum(all_dat.iloc[6:21,0]) # doublecheck
    f2 = (all_dat.loc['ABBA'] - all_dat.loc['ABAB']) / total_inv
    f1 = (all_dat.loc['AABB'] - all_dat.loc['ABAB']) / total_inv
    
    muf2 = 1 # reproduced from https://github.com/lkubatko/HilsTest/blob/master/SummaryScripts/HybTest.R
    muf1 = 0
    
    sif12 = 1/total_inv * ((all_dat.loc['AABB']/total_inv) * \
                           (1 - all_dat.loc['AABB']/total_inv) + \
                           (all_dat.loc['ABAB']/total_inv) * \
                           (1 - all_dat.loc['ABAB']/total_inv) + \
                           2 * (all_dat.loc['AABB']/total_inv) * \
                           (all_dat.loc['ABAB']/total_inv))
    sif22 = 1/total_inv * ((all_dat.loc['ABBA']/total_inv) * \
                           (1 - all_dat.loc['ABBA']/total_inv) + \
                           (all_dat.loc['ABAB']/total_inv) * \
                           (1 - all_dat.loc['ABAB']/total_inv) + \
                           2 * (all_dat.loc['ABBA']/total_inv) * \
                           (all_dat.loc['ABAB']/total_inv))
    sif1f2 = 1/total_inv * (-(all_dat.loc['AABB']/total_inv) * \
                           (all_dat.loc['ABBA']/total_inv) + \
                           (all_dat.loc['AABB']/total_inv) * \
                           (all_dat.loc['ABAB']/total_inv) + \
                           (all_dat.loc['ABBA']/total_inv) * \
                           (all_dat.loc['ABAB']/total_inv) + \
                           (all_dat.loc['ABAB']/total_inv) * \
                           (1 - all_dat.loc['ABAB']/total_inv))
        
    hils = f2*(f1/f2 - muf1/muf2) / (sif22 * (f1/f2)**2 - 2*sif1f2*(f1/f2) + sif12)**.5
    
    return hils


# changing the labels of te output figures
def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('hybrid taxon')

# make a violin-plot for a given dataframe
def play_the_violin(all_dat,
                    m_dat,
                    mapfile):
    from matplotlib import pyplot as plt

    ANALYSIS = [r'Hybridization index $\gamma$', 'D statistic', 'Hils statistic']

    fig, axes = plt.subplots(1, 3, figsize=(18, 7))
    both = list(all_dat.loc['P1'].values) + list(all_dat.loc['P2'].values)
    pars = list(sorted(set(both), key=both.index))
    hybs = list(set(all_dat.loc['Hybrid_pop']))
    paroverall = [sum(m_dat.loc['pop'] == sp) for sp in pars]
    hylength = [sum(all_dat.loc['Hybrid_pop'] == sp) for sp in hybs]
    hyoverall = [sum(m_dat.loc['pop'] == sp) for sp in hybs]
    labels = [f'{hb}\nN = {hl} ({ho})' for hb, hl, ho in zip(hybs, hylength, hyoverall)]
    
    for co, test in enumerate(['Gamma', 'D_stat', 'Hils']):
        data = []
        for hyb in hybs:
            data.append([float(x) for x in all_dat.loc[test, all_dat.loc['Hybrid_pop'] == hyb].values])
        #        [float(x) for x in all_dat.loc['Gamma', all_dat.loc['Hybrid_pop'] == 'vvwest'].values]]
        r = axes[co].violinplot(dataset = data)
        
        COLOR = ['blue', 'red', 'green']
        
        for pc in r['bodies']:
            pc.set_facecolor(COLOR[co])
        r['cbars'].set_color(COLOR[co])
        r['cmaxes'].set_color(COLOR[co])
        r['cmins'].set_color(COLOR[co])
    
        ylabs = [f'Hybr. index $\gamma$ in % genomic contribution of {pars[0]}', 
                 f'D statistic (p1: {pars[0]}, p2: {pars[1]})', 
                 f'Hils statistic (p1: {pars[0]}, p2: {pars[1]})']
        axes[co].set_title(ANALYSIS[co])
        axes[co].yaxis.grid(True)
        axes[co].set_ylabel(ylabs[co])
        set_axis_style(axes[co], labels)
    
    axes[0].set_ylim(0, 1)
    axes[1].set_ylim(0, 1)
    
    fig.suptitle(f'Theoretical hybridization between {pars[0]} (N = {paroverall[0]}) \
and {pars[1]} (N = {paroverall[1]})')
    fig_name = f'{MAPFILE.split("-map.txt")[0]}_{pars[0]}_{pars[1]}_hyb.pdf'
    plt.savefig(fig_name, bbox_inches='tight')

    return

# add new colums to hybridlist, which summarizes the hypothetical hybridizations
def hylicisvo(hybrid_list,
              sp12_dat,
              p1,
              p2):
    for hybrid in set(sp12_dat.loc['Hybrid_pop']):
        hyb_dat = sp12_dat.loc[:, sp12_dat.loc['Hybrid_pop'] == hybrid]
        hybrid_list.loc[f'{hybridlist.shape[0]}'] = [p1,
                                           p2,
                                           hybrid,
                                           hyb_dat.loc['Hybrid'].values.tolist(),
                                           np.mean(hyb_dat.loc['Gamma']),
                                           np.mean(hyb_dat.loc['D_stat']),
                                           np.mean(hyb_dat.loc['Hils'])]
    return hybrid_list


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="infile", required=True,
                    help="output-file of a HyDe run.")
parser.add_argument("-m", action="store", dest="mapfile", required=True,
                    help="map-file which was used for the HyDe run.")
args = parser.parse_args()

INFILE = args.infile
MAPFILE = args.mapfile
# INFILE = '2020_04_vitis_imputed-ind.txt'
# MAPFILE = '2020_04_vitis_imputed-map.txt'
HYB_COLUMS = ['p1', 'p2', 'hyp_hybrid', 'individuals', 'mean_gamma', 'mean_d', 'mean_hils']

d_dat = to_pededeff(INFILE, head=True)

d_dat[3:6] = d_dat.iloc[3:6].astype('float')
d_dat[6:] = d_dat.iloc[6:].astype('float')
d_dat[6:] = d_dat.iloc[6:].astype('int')

d_dat.loc['D_stat'] = mama_mia(d_dat)
d_dat.loc['Hils'] = up_the_hils(d_dat)

sign = [x for x, p in enumerate(d_dat.loc['Pvalue']) if p <= .05]
d_sign = d_dat.iloc[:,sign]

specs = set(pd.concat([d_sign.loc['P1'],d_sign.loc['P2']]))

all_dat = d_sign.loc[:, list(d_sign.loc['P1'] == 'nix')].copy()


for spo, spe in combinations(specs, 2):
    subs = d_sign.loc[:, list(d_sign.loc['P1'] == spo)]
    sabs = d_sign.copy().loc[:, list(d_sign.loc['P2']==spo)]
    # sabs.loc['Gamma'] = 1 - sabs.loc['Gamma'] # 1 - lambda!!! - not approproiate! would need change of site patterns!
    subss = subs.loc[:, list(subs.loc['P2'] == spe)]
    sabss = sabs.loc[:, list(sabs.loc['P1'] == spe)] 
    sebss = sabss.copy() # change all P1 to same taxon and all P2 to same taxon. -> doubled values in df!
    # sebss.loc['P2'] = sabss.copy().loc['P1']
    # sebss.loc['P1'] = sabss.copy().loc['P2']
    if subss.values.size != 0:
        all_dat = pd.concat([all_dat, subss], axis=1)
    elif sabss.values.size != 0:
        all_dat = pd.concat([all_dat, sabss], axis=1)
    # sebbo = pd.concat([subss,sebss], axis=1)
    # all_dat = pd.concat([all_dat, sebbo], axis=1)


m_dat = to_pededeff(MAPFILE, head=False)
spacs = set(m_dat.loc['pop'])

hybrid_line = pd.Series([""] * all_dat.shape[1])
all_dat.loc['Hybrid_pop'] = hybrid_line

for spo in spacs:
    spo_ind = m_dat.loc['ind'][m_dat.loc['pop'] == spo]
    plot_sep = all_dat.loc['Hybrid'].isin(spo_ind)
    plot_frame = all_dat.loc[:, plot_sep]
    all_dat.loc['Hybrid_pop', plot_sep] = [spo] * sum(plot_sep)

hybridlist = pd.DataFrame(columns = HYB_COLUMS)
for sp1, sp2 in combinations(specs, 2):
    sp1_dat =  all_dat.loc[:,[v for v in all_dat.loc['P1'] == sp1]]
    sp12_dat = sp1_dat.loc[:,[v for v in sp1_dat.loc['P2'] == sp2]]
    if sp12_dat.values.size != 0:
        play_the_violin(sp12_dat, m_dat, MAPFILE)
        hybridlist = hylicisvo(hybridlist, sp12_dat, sp1, sp2)
    else:
        sp2_dat =  all_dat.loc[:,[v for v in all_dat.loc['P1'] == sp2]]
        sp12_dat = sp2_dat.loc[:,[v for v in sp2_dat.loc['P2'] == sp1]]
        if sp12_dat.values.size != 0:
            play_the_violin(sp12_dat, m_dat, MAPFILE)
            hybridlist = hylicisvo(hybridlist, sp12_dat, sp2, sp1)

hybridlist.to_csv(f'{MAPFILE.split("-map.txt")[0]}_hyde_hybr.csv')
