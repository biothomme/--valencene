#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 08:25:18 2020

this script is thought to permutate all individual of 5 populations with each other.
that can be useful for Dfoil analysis...

@author: Thomsn
"""
import argparse
import itertools as it
import os
import pandas as pd

def all_ind_per_tax(indiv_dat, taxa):
    get_low = lambda x: ''.join(''.join(x.lower().split('.')).split(' '))
    indiv_dat['Pop-celine'] = [get_low(ta) for ta in indiv_dat['Pop-celine']]
    tax_ass = pd.Series(index = taxa)
    for tax in taxa:
        ind_list = []
        if tax != 'vog':
            for _, ind in indiv_dat.iterrows():
                if ind['Pop-celine'] == tax:
                    ind_list.append(ind['cut_names'])
        else:
            for _, ind in indiv_dat.iterrows():
                if ind['Pop-celine'] not in taxa[:4]:
                    ind_list.append(ind['cut_names'])
        tax_ass[tax] = ind_list
    return tax_ass

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="summary_file", required=True,
                    help="Name of the input fasta file")


args = parser.parse_args()
IND_SUMMARY = args.summary_file

# WORDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
# IND_SUMMARY = '2020_04_vitis_USA_2_summary.csv'
TAXA = ['vveast',
        'vsylveast',
        'vvwest',
        'vsylvwest',
        'vog']

# os.chdir(WORDIR)
ind_dat = pd.read_csv(IND_SUMMARY)
taxon_assignments = all_ind_per_tax(ind_dat, TAXA)

output_file = f'{"/".join(IND_SUMMARY.split("/")[:-1])}/2020_04_vitis_5tax_perm.csv'
with open(output_file, 'w') as of:
    of.write(f'{",".join(TAXA)}\n')
    for coco, comb in enumerate(it.product(taxon_assignments[0],
                           taxon_assignments[1],
                           taxon_assignments[2],
                           taxon_assignments[3],
                           taxon_assignments[4])):
        of.write(f'{",".join(comb)}\n')
        print(coco)
