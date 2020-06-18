#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 14:06:03 2020

script to obtain the best networks from a phylonet
output.

@author: Thomsn
"""

import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="bigdir", required=True,
                    help="bigdir with repdirs")

args = parser.parse_args()
WORDIR = args.bigdir
TREEDIR = 'seed_rep_new_mnr1_ret0'

dirro = [elm for elm in os.listdir(WORDIR) if 'rep' in elm]
back_dirs = len(dirro[0].split('/')) + len(WORDIR.split('/'))

back_door = '/'.join(['..']*back_dirs)

for dirre in dirro:
    # print(dirre)
    topoframe = pd.DataFrame(columns = ['pseudo_likelihood', 'topology'])
    cont = [repl for repl in os.listdir(f'{WORDIR}/{dirre}/{TREEDIR}') if 'net0.nex' in repl]
    for topfile in cont:
        with open(f'{WORDIR}/{dirre}/{TREEDIR}/{topfile}', 'r') as tfl:
            soup = tfl.readlines()
            now = False
            for molecule in soup:
                if now:
                    loglik = molecule.split(':')[0]
                    topomol = molecule.split(']')[1]
                    now = False
                    seed_num = topfile.split('_ret_')[1].split('.')[0]
                    topoframe.loc[seed_num] = [loglik, topomol]
                elif 'Topology' in molecule:
                    now = True
    topoframe = topoframe.sort_values(by=["pseudo_likelihood"])
    topoframe.to_csv(f'{WORDIR}/{dirre}/tree_sum.csv')



