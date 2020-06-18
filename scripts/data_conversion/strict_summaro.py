#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 07:47:30 2020

this script may be used to obtain a strict summary of unbiased
 ind of the vitis dataset

@author: Thomsn
"""

import os
import pandas as pd

WORDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/additional_information/structure'
INFILE_SYLV = 'Structure-SYLVESTRIS-v1.csv'
INFILE_SAT = 'Structure-SATIVA-v1.csv'
IND_SUMMARY = '../../2020_04_vitis_USA_2_summary.csv'
SUSP_IND = '../../hyde/admixed_asia/suspiciousindividuals.csv'
SYLTL = {'wceu': 'vsylvwest',
         'emca': 'vsylveast'}
VINTL = {'wceu': 'vvmeds',
         'mfea': 'vveast',
         'balk': 'vvbalk'}

os.chdir(WORDIR)

ind_dat = pd.read_csv(IND_SUMMARY)
susp_dat = pd.read_csv(SUSP_IND, sep=';')

strict_ind = []
strict_tax = []
for spc, species in enumerate([INFILE_SYLV, INFILE_SAT]):
    ## first assign any individual to a group 
    struct_dat = pd.read_csv(species)
    if spc == 0:
        TAXA = ['vsylveast', 'vsylvwest']
        struct_dat = struct_dat.copy()[struct_dat['ADN-ID'] != 'Vbg152']
    else:
        TAXA = ['vvwest', 'vveast']
    
    for _, indi in struct_dat.iterrows():
        if spc == 1:
            # criteria = pd.Series(data = [indi.loc['Structure-A1 BALK']+indi.loc['Structure-A2 WCEU'],
            #       indi.loc['Structure-A3 MFEA']],
            #       index = ['balk+wceu',
            #                'mfea'])
            criteria = pd.Series(data = [indi.loc['Structure-A1 BALK'],
                                         indi.loc['Structure-A2 WCEU'],
                                         indi.loc['Structure-A3 MFEA']],
                  index = ['balk',
                           'wceu',
                           'mfea'])
            sig_grp = 'mfea'
        else:
            criteria = pd.Series(data = [indi.loc['Structure A1 EMCA'],
                  indi.loc['Structure A2 WCEU']],
                  index = ['emca',
                           'wceu'])
            sig_grp = 'wceu'
        save_drop = criteria[criteria > .85]
        if len(save_drop) > 0:
            # print(criteria)
            nm = save_drop.index[0]
            if spc == 1:
                taxo = VINTL[nm]
            else:
                taxo = SYLTL[nm]
            strict_tax.append(taxo)
            strict_ind.append(indi['ADN-ID'])

strict_dat = pd.DataFrame(columns=list(ind_dat.columns)+['taxnew'])
for _, indi in ind_dat.iterrows():
    if (indi['ADN-ID'] in strict_ind):
        k = [i for i, val in enumerate(strict_ind) if val == indi['ADN-ID']][0]
        strict_dat.loc[str(len(strict_dat))] = list(indi.values) + [strict_tax[k]]
    elif (indi['ADN-ID'] not in list(susp_dat['ADN-ID'].values)) & \
    (indi['Pop-celine'] in ['Vitis Asia', 'Vitis USA', 'Vitis USA2']):
        strict_dat.loc[str(len(strict_dat))] = list(indi.values) + [indi['Pop-celine'].lower().replace(' ','')]

strict_dat.to_csv('../../2020_04_vitis_strc_smry_balk.csv')
 
 