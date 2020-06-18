#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 11:44:59 2020

this seems to be the seed of a crazy project. this script should enable to
split robs excel-file into different files, which are needed for KDE analysis.

@author: Thomsn
"""
import os
import pandas as pd
import random

WORDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
EXCEL_FILE = '2020_04_vitis_USA_2.txt'
GEN_SUM = '2020_04_vitis_USA_2_summary.csv'
STRC_SUM = '2020_04_vitis_strc_smry_balk.csv'

os.chdir(WORDIR)

shortify = lambda x: ''.join(x.split('_'))[-10:].lower()
SNP_SET = 'snp_set.csv'
HETEROMAX = .5

with open(EXCEL_FILE, 'r') as ef:
#    with open(SNP_SET, 'w') as ss:
        bim_list = pd.DataFrame(columns = ['chr', 'pnt', 'cm', 'bp', 'ref', 'alt'])
        soup = ef.readlines()
        for co, molecule in enumerate(soup):
            if co == 0:
                all_ind = [shortify(ii) for ii in molecule.split('\n')[0].split('\t')[1:]]
                # ss.write(','.join(all_ind))
                # ss.write('\n')
            elif co > 4:
                line = molecule.split('\n')[0]
                chrm = line.split('chr')[1].split('_')[0]
                crd = line.split('_')[1].split('\t')[0]
                all_bases = ''.join(line.split('\t')[1:])
                variants = list(set(all_bases))
                if len(variants) > 2:
                    break
                if chrm == 'Un':
                    chrm = '99'
                if crd == 'random':
                    crd = random.randint(10000, 1000000)
                    chrm = '99'
                dominance = [all_bases.count(vari) for vari in variants]
                ref_alt_var = [va for _, va in sorted(zip(dominance,variants), reverse = True)]
                bim = [int(chrm),
                       '.',
                       '0',
                       int(crd),
                       ref_alt_var[0],
                       ref_alt_var[1]]
                bim_list.loc[f'{chrm}_{crd}'] = bim
                # bp_list = []
                # for coco, basepair in enumerate(line.split('\t')[1:]):
                #     indici = list(snp_frame.copy().index) + [f'{chrm}_{crd}']
                #     bp_list.append(str(basepair.count(ref_alt_var[1])))
                # ss.write(','.join(bp_list))
                # ss.write('\n')
snp_set = pd.read_csv(SNP_SET)
snp_set.index = bim_list.index
bim_list = bim_list.sort_values(by=['chr','bp']).copy()
# new_snp_frame = pd.DataFrame(columns = snp_frame.columns)
# for _, bim in bim_list.iterrows():
#     try:
#         new_snp_frame.loc[bim.name] = snp_set.loc[bim.name,:]
#     except KeyError:
#         pass
# new_snp_frame.to_csv('sorted_snp_set.txt')
new_snp_frame = pd.read_csv('sorted_snp_set.txt')
new_snp_frame.index = bim_list.index
new_snp_frame = new_snp_frame.drop(columns='Unnamed: 0').copy()
nsf_len = new_snp_frame.shape
heterozygosity = [sum([i for i in fr if i == 1])/nsf_len[1] for _, fr in new_snp_frame.iterrows()]
heterozygosity = pd.Series(data=heterozygosity,
                              index=new_snp_frame.index)
indheterozygosity = [sum([i for i in new_snp_frame[fr] if i == 1])/nsf_len[0] for fr in list(new_snp_frame)]
indheterozygosity = pd.Series(data=indheterozygosity,
                              index=new_snp_frame.columns)
indheterozygosity[indheterozygosity < HETEROMAX]
selected_snps = new_snp_frame[heterozygosity < HETEROMAX].copy()
bim_selected = bim_list[heterozygosity < HETEROMAX].copy()

with open('vitis.bim', 'w') as vb:
    for _, line in bim_selected.iterrows():
        vb.write('\t'.join([str(x) for x in line.values]))
        vb.write('\n')
        
with open('vitis.fam', 'w') as vf:
    for i, ind in enumerate(all_ind):
        fam = [ind, ind, '0', '0', '0', '-9']
        vf.write('\t'.join(fam))
        vf.write('\n')

with open('vitis.geno', 'w') as vg:
    for col in selected_snps.columns:
        vg.write(''.join([str(x) for x in selected_snps[col].values]))
        vg.write('\n')

gen_sum = pd.read_csv(GEN_SUM).drop(columns = ['Unnamed: 0'])
strc_sum = pd.read_csv(STRC_SUM).drop(columns = ['Unnamed: 0', 'Unnamed: 0.1'])
gen_sum['Pop-celine'] = [x.replace('.','').replace(' ','').lower() for x in gen_sum['Pop-celine']]

taxid = 1
with open('vitis.ref', 'w') as vr:
    with open('vitis_taxon_id.txt', 'w') as ti:
        for taxx in set(strc_sum['taxnew']):
            ti.write(f'{taxid}\t{taxx}\n')
            sub_sum = strc_sum[strc_sum['taxnew'] == taxx]
            for _, line in sub_sum.iterrows():
                vr.write(f'{taxid}\t{line["cut_names"]}')
                vr.write('\n')
            taxid += 1

with open('vitis.admix', 'w') as va:
    with open('vitis_taxon_id.txt', 'a') as ti:
        for taxx in set(gen_sum['Pop-celine']):
            ti.write(f'{taxid}\t{taxx}_admix\n')
            sub_sum = gen_sum[gen_sum['Pop-celine'] == taxx]
            for _, line in sub_sum.iterrows():
                if line["cut_names"] not in list(strc_sum["cut_names"].values):
                    va.write(f'{taxid}\t{line["cut_names"]}')
                    va.write('\n')
            taxid += 1