#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 07:17:56 2020

this script is thought to open phylonet nexus files and randomly change 
the seeds. output are seed replicate files.

@author: Thomsn
"""
import argparse
import os


# functions
def new_nex(infile, topology, outfile):
    import random
    MNR = 1
    MEC = 500000
    MF = 500
    PL = 10
    MR = 1
    seed = random.randint(1, 10000000)
    with open(topology, 'r') as topol:
        soup = topol.readlines()
        for i, molecule in enumerate(soup):
            if i == 0:
                topo_str = molecule
    with open(infile, 'r') as infl:
        with open(outfile, 'a') as oufl:
            soup = infl.readlines()
            for molecules in soup:
                if '-pseudo' in molecules:
                    tax_str = molecules.split('-tm ')[1].split('>')[0]
                    seedstr = f'MLE_BiMarkers -pseudo -mnr {MNR} -mec {MEC} \
-mf {MF} -sd {seed} -pl {PL} -mr {MR} -tm {tax_str}> -esptheta -snet (net1);\n'
                    oufl.write(seedstr)
                elif 'Network net1 = ' in molecules:
                    seedstr = f'Network net1 = {topo_str}\n'
                    oufl.write(seedstr)
                else:
                    oufl.write(molecules)
    print(f'{outfile} with new seed was produced...')
    return 


# WORDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonet/2020_04_vitis_USA_2_summary_20_reps/rep_01'
# NEXFILE = 'rep_std_data_0_ret.nexus'

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', 
                    dest = 'nexfile',
                    action='store', 
                    required=True,
                    help='this needs to be a phylonet nexus file with the word -pseudo in the phylonet command.')
parser.add_argument('-t',
                    dest = 'newick',
                    action='store', 
                    required=True,
                    help='this needs to be the new topology file.')

args = parser.parse_args()
NEXFILE = args.nexfile
NEWICK = args.newick
MR = 1
# os.chdir(WORDIR)
# stem = NEXFILE.split('.nexus')[0].split('/')[-1]
stem = f"_{MR}_".join(NEXFILE.split(".nexus")[0].split("_0_")).split('/')[-1]
n_files_exist = sum([1 for fl in os.listdir('/'.join(NEXFILE.split('/')[0:-1])) if stem in fl])

outnex = f'{f"_{MR}_".join(NEXFILE.split(".nexus")[0].split("_0_"))}_{n_files_exist}.nexus'

new_nex(NEXFILE, NEWICK, outnex)
