#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 09:00:59 2020

This script is thought to produce replicates of SNP datasets (here: vitis)
by using structure analysis probabilities as well as ibs clusters to avoid
a close realtionships

@author: Thomsn
"""
# choices for sativa
# dataset: vveast = A3 >= .75
# dataset: vvwest = A3 < .75

__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'

import os
import numpy as np
import pandas as pd
import random
import argparse

TOTAL_REP = 20
NUM_EUR = 5 # how many samples each of vveast, vvwest, vsylveast, vslylvwest
NUM_ASIA = 1
FRESHOLD = .2
VITISASIA = ['V. piasezkii',
             'V. armata',
             'V. davidii',
             'V. flexuosa',
             'V. romanetii',
             'V. balansaeana']

VITISONEASIA = ['V. yeshanensis',
                'V. coignetiae',
                'V. thunbergii',
                'V. pentagona']

VITISUSA = ['V. muscadinia']


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-y", action="store", dest="infile_sylv", required=True,
                    help="structure excel sheet for sylvestris")
parser.add_argument("-a", action="store", dest="infile_sat", required=True,
                    help="structure excel sheet for sativa")
parser.add_argument("-c", action="store", dest="clusters", required=True,
                    help="file containing all assigned clsuters as csv")
parser.add_argument("-i", action="store", dest="ind_summary", required=True,
                    help="summary csv file with all individuals used")

args = parser.parse_args()
infile_sylv = args.infile_sylv
infile_sat = args.infile_sat
clusters = args.clusters
ind_summary = args.ind_summary

# wordir = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/additional_information/structure'
# infile_sylv = 'Structure-SYLVESTRIS-v1.csv'
# infile_sat = 'Structure-SATIVA-v1.csv'
# clusters = '../ibs_cluster/clusters_vitis_13sativa_10sylvestris.csv'
# ind_summary = '../../2020_04_vitis_USA_summary.csv'

# os.chdir(wordir)

ind_dat = pd.read_csv(ind_summary)
clus_dat = pd.read_csv(clusters)
st1 = ['_'.join(grob.split('GT')[0].split('x')) for grob in clus_dat['adn_id']] # reconvert the x_names
clus_dat['adn_id'] = [f'{grob}GT' if 'Guruchay_' in grob else grob for grob in st1]

rep_path = f'{ind_summary.split(".csv")[0]}_{TOTAL_REP}_reps'
if not os.path.exists(rep_path):
    os.makedirs(rep_path)



## get replicates for the species v. sylvestris and v. sativa
for rep_num in range(TOTAL_REP):
    repl_list = []
    tclust_list = []
    if len(str(rep_num + 1)) == 1:
        nr = f'0{rep_num + 1}'
    else:
        nr = f'{rep_num + 1}'
    new_path = f'{rep_path}/rep_{nr}'
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    for species in [infile_sylv, infile_sat]:
        struc_dat = pd.read_csv(species)
        #x_the_uscr = ['x'.join(adname.split('GT')[0].split('_')) for adname in struc_dat['ADN-ID'].values]
        for taxon in ['west', 'east']:
            clus_list = []
            sm_repl_list = []
            if species == infile_sat:
                if taxon == 'east':
                    slc_struc_dat = struc_dat[struc_dat['Structure-A3 MFEA'] > FRESHOLD]
                    rand_frag = slc_struc_dat['Structure-A3 MFEA'] / sum(slc_struc_dat['Structure-A3 MFEA'])
                else:
                    slc_struc_dat = struc_dat[struc_dat['Structure-A3 MFEA'] <= (1 - FRESHOLD)]
                    rand_frag = (1 - slc_struc_dat['Structure-A3 MFEA']) / sum(1 - slc_struc_dat['Structure-A3 MFEA'])
            else:
                if taxon == 'east':
                    slc_struc_dat = struc_dat[struc_dat['Structure A1 EMCA'] > FRESHOLD]
                    rand_frag = slc_struc_dat['Structure A1 EMCA'] / sum(slc_struc_dat['Structure A1 EMCA'])
                else:
                    slc_struc_dat = struc_dat[struc_dat['Structure A1 EMCA'] <= (1 - FRESHOLD)]
                    rand_frag = (1 - slc_struc_dat['Structure A1 EMCA']) / sum(1 - slc_struc_dat['Structure A1 EMCA'])            
            rand_template = np.cumsum(rand_frag)
            
            roun = 1
            while len(sm_repl_list) < NUM_EUR:
                rnum = random.random()
                choice = int([i for i, value in enumerate(rand_template) if rnum < value][0])
                selected_ind = slc_struc_dat.iloc[choice,:].loc['ADN-ID']
                selected_cluster = clus_dat[clus_dat['adn_id'] == selected_ind]['ward.D2'].values[0]
                if selected_ind not in sm_repl_list:
                    if selected_ind not in repl_list:
                        if sum([True for cl in clus_list if cl == selected_cluster]) < int(np.log10(roun) / np.log10(100) -1):  # to not get stuck in the loop
                            sm_repl_list.append(selected_ind)
                            clus_list.append(selected_cluster)
                roun += 1
            repl_list += sm_repl_list
            tclust_list += clus_list

    ## get replicates for the asian taxa
    ## 1 ind each of V. piasezkii, V. armata, V. davidii, V. flexuosa, V. romanetii, V. balansaeana
    ## 1 ind of V. muscadinia
    ## 1 ind of the group V. yeshanensis, V. coignetiae, V. thunbergii, V. pentagona
    
    ind_dat['taxon'] = [' '.join(tax.split(' ')[0:2]) for tax in ind_dat['taxon']]  #get rid of blank spaces after taxon names
    ind_dat['taxon'] = ['V. balansaeana' if tax == 'V. balanseana' else tax for tax in ind_dat['taxon']] # get rid of two balansaeanas with different spellings
    set(ind_dat['taxon'])
    
    
    for spec in VITISASIA + VITISUSA + [VITISONEASIA]:
        if len(spec) == 1:
            sep_dat = ind_dat[ind_dat['taxon'] == spec]
        else:
            sep_bool = [True if sp in spec else False for sp in ind_dat['taxon']]
            sep_dat = ind_dat[sep_bool]
        offer = sep_dat['ADN-ID'].values
        pick = random.sample(list(offer), NUM_ASIA)
        repl_list += pick
        tclust_list += ['NA']
    
    ## subset the datasets for the sampled ind
    repl_bool = [True if adn in repl_list else False for adn in ind_dat['ADN-ID']]
    repl_dat = ind_dat.copy()[repl_bool]
    clustri = pd.Series(data = tclust_list, index = repl_list)
    repl_dat['cluster'] = clustri[repl_dat['ADN-ID']].values
    repl_dat = repl_dat.drop('Unnamed: 0',1)
    sorted(repl_list)
    repl_dat.to_csv(f'{new_path}/replicate_data.csv')
