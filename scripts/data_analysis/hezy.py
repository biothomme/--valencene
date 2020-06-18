#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:27:58 2020

@author: Thomsn
"""
from matplotlib import pyplot as plt
import numpy as np
import os
import pandas as pd

BIGDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
SNPS = 'snp_set.csv'
SUMRY = '2020_04_vitis_USA_2_summary.csv'

os.chdir(BIGDIR)
os.listdir()

snps = pd.read_csv(SNPS)
smry = pd.read_csv(SUMRY)

plot_frame = pd.DataFrame(columns=['taxon', 'heterozygosity', 'name'])
for tax in sorted(set(smry['Pop-celine'])):
    tox = tax
    for prfx in ['ia', 'itis', ' ', '.', 'est', 'ast', 'ylv']:
        tax = tax.replace(prfx, '')
    tax = tax.lower().replace('usa', 'us').replace('us2', 'u2')
    print(tax)
    taxset = smry.loc[smry['Pop-celine'] == tox,:]
    hezy = []
    inds = []
    for ind in taxset['cut_names']:
        hezy.append(sum([1 for snp in snps[ind] if snp == 1])/len(snps[ind]))
        inds.append(ind)
    plot_frame.loc[len(plot_frame),:] = [tax, hezy, inds]


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Taxonomic group')

# make a violin-plot for a given dataframe


ANALYSIS = [r'Hybridization index $\gamma$', 'D statistic', 'Hils statistic']

fig, axes = plt.subplots(1, 1, figsize=(8, 5))
r = axes.violinplot(dataset = plot_frame['heterozygosity'])
COLOR = ['blue', 'red', 'green']

for pc in r['bodies']:
    pc.set_facecolor(COLOR[2])
r['cbars'].set_color(COLOR[2])
r['cmaxes'].set_color(COLOR[2])
r['cmins'].set_color(COLOR[2])
ylabs = 'Percentage of heterozygosity'
axes.yaxis.grid(True)
axes.set_ylabel(ylabs)
set_axis_style(axes, plot_frame['taxon'])


fig_name = f'heterozgosity.pdf'
plt.savefig(fig_name, bbox_inches='tight')
