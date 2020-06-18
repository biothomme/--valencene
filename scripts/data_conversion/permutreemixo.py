#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:25:15 2020

this script should convert a biallelic matrix and population 
assignment replicates (from ibs clusters) to get an output that can be fed into treemix

@author: Thomsn
"""

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="big_dir", required=True,
                    help="Name of the input big dir")
parser.add_argument('-n', action='store', dest='rep_dir',
                    help='dir with replicates inside')

args = parser.parse_args()
BIG_DIR = args.big_dir
REP_DIR = args.rep_dir

# BIG_DIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
# REP_DIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/treemix/permutations'
SNP_FILE = 'sorted_snp_set.txt'
POP_FILE = '2020_04_vitis_USA_2_summary.csv'
SUSP_IND = ['b00ei5i',
             'b00er6s',
             'vbg267',
             'b00erc1',
             'b00erco',
             'vbg90',
             'b001757',
             'vbg180',
             'vbg88',
             'b00fqik',
             'vbg175',
             'vbg152']

rep_post = '/'.join([path for path in REP_DIR.split('/') if path not in BIG_DIR.split('/')])
# os.chdir(BIG_DIR)
rep_list = [rep for rep in os.listdir(f'{BIG_DIR}/{rep_post}') if 'vs_' in rep]

snp_array = pd.read_csv(f'{BIG_DIR}/{SNP_FILE}')
snp_array.index = list(snp_array.iloc[:,0].values)
snp_array = snp_array.drop(columns='Unnamed: 0')

pop_array = pd.read_csv(f'{BIG_DIR}/{POP_FILE}')
pop_array.index = list(pop_array.iloc[:,0].values)
pop_array = pop_array.drop(columns='Unnamed: 0')
cutnms = lambda x: x.lower().replace(' ','').replace('.','')
pop_array['Pop-celine'] = [cutnms(x) for x in pop_array['Pop-celine']]


# os.chdir(rep_post)
for repp in rep_list:
    # os.chdir(repp)
    cluster_file = [fi for fi in os.listdir(f'{BIG_DIR}/{rep_post}/{repp}') if '_vitis' in fi][0]
    clu_array = pd.read_csv(f'{BIG_DIR}/{rep_post}/{repp}/{cluster_file}')
    std_names = [x.replace('x','')[-8:].lower() for x in clu_array['adn_id']]
    std_names = [f'{x}gt' if 'ay' in x else x for x in std_names]
    clu_array['adn_id'] = std_names
    clu_dat = pd.DataFrame(columns=['tax', 'regions'])
    for _, row in clu_array.iterrows():
        cl = row['ward.D2']
        if cl not in list(clu_dat.index):
            clu_dat.loc[cl,:] = [row['spnames'],[row['region']]]
        else:
            if row['region'] not in clu_dat.loc[cl,'regions']:
                clu_dat.loc[cl,'regions'] = clu_dat.loc[cl,'regions'] + [row['region']]
    maxname = max([len(str(x)) for x in clu_dat.index])
    clu_names = []
    for _, cd in clu_dat.iterrows():
        if cd['tax'] == 'sativa':
            tx = 'VV'
        else:
            tx = 'VS'
        cdreg = cd['regions']
        cdreg.sort()
        regs = '_'.join([rg[:2].lower() for rg in cdreg])
        cdn = ''.join((maxname - len(str(cd.name))) * ['0'] + [str(cd.name)])
        clu_names.append(f'{cdn}_{tx}_{regs}')
    clu_dat['names'] = clu_names
    
    clu_names = []
    cla_names = []
    for _, row in pop_array.iterrows():
        if row['cut_names'] in list(clu_array['adn_id'].values):
            curar = clu_array.loc[clu_array['adn_id'] == row['cut_names'],:]
            clus = curar['ward.D2'].values[0]
            clus_name = clu_dat.loc[clus,'names']
            clas_name = clus_name
        else:
            clus_name = pop_array.loc[pop_array['cut_names'] == row['cut_names'],:]['taxon'].values[0]
            clus_name = clus_name.lower().replace('.', '').replace(' ','')
            clas_name = pop_array.loc[pop_array['cut_names'] == row['cut_names'],:]['taxon'].values[0]
            clas_name = clas_name.lower().replace('.', '').replace(' ','')
        clu_names.append(clus_name)
        cla_names.append(clas_name)
    pop_array['cluster_og'] = clu_names
    pop_array['cluster'] = cla_names

    with open(f'{BIG_DIR}/{rep_post}/{repp}/vitis_treemix_og.txt', 'w') as vt:
        header = []
        for pop in set(pop_array['cluster_og'].values):
            header.append(pop)
        vt.write(f'{" ".join(header)}\n')
        print(repp)
        for _, line in snp_array.iterrows():
            pop_counts = []
            for pop in set(pop_array['cluster_og'].values):
                inds = list(pop_array[pop_array['cluster_og'].values == pop]['cut_names'].values)
                inds = [ind for ind in inds if ind not in SUSP_IND]
                total_inds = len(inds) * 2
                alleles = line[inds]
                nonal = sum(alleles)
                domal = total_inds - nonal
                pop_counts.append(f'{domal},{nonal}')
            vt.write(f'{" ".join(pop_counts)}\n')
            
    with open(f'{BIG_DIR}/{rep_post}/{repp}/vitis_treemix.txt', 'w') as vt:
        header = []
        for pop in set(pop_array['cluster'].values):
            header.append(pop)
        vt.write(f'{" ".join(header)}\n')
            
        for _, line in snp_array.iterrows():
            pop_counts = []
            for pop in set(pop_array['cluster'].values):
                inds = list(pop_array[pop_array['cluster'].values == pop]['cut_names'].values)
                inds = [ind for ind in inds if ind not in SUSP_IND]
                total_inds = len(inds) * 2
                alleles = line[inds]
                nonal = sum(alleles)
                domal = total_inds - nonal
                pop_counts.append(f'{domal},{nonal}')
            vt.write(f'{" ".join(pop_counts)}\n')
