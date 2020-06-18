#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 16:54:31 2020

this script is thought to force monophyly in topologies of non monophyletic taxa.
it uses output of monophorce.r and a given fasta file and produces a phylip
sequence file as well as a newick topology, which can be used in FastME for 
branch length calculations. 

input:  - bigdir with repdirs with monophorce.r output and alignments (fasta)
        - topology of a scenario which should be forced (taxon-names must fit 
          with taxa of monophorce output)

@author: Thomsn
"""
__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'

import argparse
import os
import pandas as pd
import re
from Bio import AlignIO

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="rep_dir", required=True,
                    help="Name of the bigdir")
parser.add_argument('-n', action='store', dest='topo',
                    help='newick string with scenario given')

args = parser.parse_args()
BIG_DIR = args.rep_dir
TOPO = args.topo

# BIG_DIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonet/2020_04_vitis_USA_2_summary_20_reps'
# TOPO = '../scen_class.new'
with open(TOPO, 'r') as file:
    soup = file.readlines()
    for molecule in soup:
        taxxer = re.findall('\w+', molecule)
        topology = molecule

os.chdir(BIG_DIR)
        
rep_list = [diro for diro in os.listdir() if 'rep_' in diro]
for repsdir in rep_list:
    for type_an in ['std', 'new']:
        tropology = topology
        string_collection = pd.Series(index = taxxer, dtype = str)
        for taxl in taxxer:
            topo_file = f'{repsdir}/topo_{type_an}_{taxl}.nex'
            with open(topo_file, 'r') as tfile:
                soup = tfile.readlines()
                for molecule in soup:
                    if '[&R] ' in molecule:
                        full_topo = molecule.split('[&R] ')[1]
                        half_topo = re.sub(r':([0-9])', '', full_topo).split(';')[0]
                        half_topo = re.sub(r'(\.)([0-9])*', '', half_topo)
            if len(half_topo.split(',')) < 2:
                half_topo = half_topo[1:-1]
            string_collection[taxl] = half_topo
            
        for tax in taxxer:
            tropology = (re.sub(tax, string_collection[tax], tropology))
        newick_topo = f'{repsdir}/forced_mp_wobl_{type_an}.newick'
        with open(newick_topo, 'a') as new_file:
            new_file.write(tropology)
        
        ## to phylip conversion:
        fas_file = f'{repsdir}/rep_{type_an}_asgnmt.fasta'
        alignment = AlignIO.read(fas_file, 'fasta')
        spac = []
        map_str = ''
        for _, seq in enumerate(alignment):
            eidi = seq.id.split('\\n')[0].split('|')
            eindi = ''.join(eidi[1].split('_'))[-9:].lower()
            eispi = ''.join(eidi[0].split('.')).lower()
            spac.append(eispi)
        
        phy_name = f'{repsdir}/rep_{type_an}_asgnmt.phy'

        with open(phy_name, 'a') as data_f:
            data_f.write(f'{len(spac)} {len(alignment[0].seq)}\n')
            
        spec = []
        ind = []
        for _, seq in enumerate(alignment):
            eidi = seq.id.split('\\n')[0].split('|')
            eindi = ''.join(eidi[1].split('_'))[-9:].lower()
            eispi = ''.join(eidi[0].split('.')).lower()
            spec.append(eispi)
            ind.append(eindi)
            n = 10
            physeq = ' '.join([str(seq.seq)[i:i+n] for i in range(0, len(str(seq.seq)), n)])
            # rewrite the fastafile for HyDe as phylip - names need to be shortened to last 10 chars
            with open(phy_name, 'a') as data_f:
                data_f.write(f'{eindi}{" "*int(10 - len(eindi))}')
                data_f.write(f'{physeq}\n')

