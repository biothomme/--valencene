#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 20:22:58 2020

@author: Thomsn
"""

import os
import pandas as pd

BIGDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonetworks/os3_data'

os.chdir(BIGDIR)

rep_list = [rep for rep in sorted(os.listdir()) if 'rep_' in rep]
for rep in rep_list:
    if 'pnw_out_net1.out' not in os.listdir(rep):
        networks = pd.DataFrame(columns = ['loglik', 'topo'])
        log_file = f'{rep}/pnw_out_net1.log'
        with open(log_file, 'r') as lf:
            soup = lf.readlines()
            ll = False
            for molecule in soup:
                if '-loglik=' in molecule:
                    loglik = float(molecule.split('-loglik=')[1].split(' ')[0])
                    ll = True
                if (');' in molecule and ll):
                    topo = molecule
                    networks.loc[len(networks),:] = [loglik, topo]
                    ll = False
        networks = networks.sort_values(by = 'loglik')
        print(networks)
        with open(f'{rep}/pnw_out_net1.out','w') as nf:
            nf.write(f'{networks["topo"][0]} -Ploglik = {networks["loglik"][0]}\n')
            nf.write('-------\nList of estimated networks for all runs (sorted by log-pseudolik; the smaller, the better):\n')
            for _, row in networks.iterrows():
                print(f' {row["topo"]} with -loglik {row["loglik"]}\n')
                nf.write(f' {row["topo"]} with -loglik {row["loglik"]}\n')
            nf.write('-------\n')
                             
                         
                         
                         
                         