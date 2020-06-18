#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:33:08 2020

this script is made to split a fasta dataset into a smaller one,
by using an input selection as filter

@author: Thomsn
"""
__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'


import os
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="infile", required=True,
                    help="total fasta file")
parser.add_argument("-d", action="store", dest="rep_dir", required=True,
                    help="path of the directory with the rep directories")

# wordir = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
# infile = '2020_04_vitis_USA.txt.fasta'
# rep_dir = 'phylonet/2020_04_vitis_USA_summary_20_reps'

args = parser.parse_args()
infile = args.infile
rep_dir = args.rep_dir

rep_name = f'{infile.split("/")[-1].split(".")[0]}_replicate.fasta'

# os.chdir(wordir)
rep_list = [diro for diro in os.listdir(rep_dir) if 'rep_' in diro]
for rep in rep_list:
    sel_file = f'{rep_dir}/{rep}/replicate_data.csv'
    fasta_rep = f'{rep_dir}/{rep}/rep_std_asgnmt.fasta'
    fasta_new = f'{rep_dir}/{rep}/rep_new_asgnmt.fasta'  # this is the file with new assignments
    alignment = AlignIO.read(infile, 'fasta')
    sel_dat = pd.read_csv(sel_file)
    for seq in alignment:
        num = [(nu, val) for nu, val in enumerate(sel_dat['ADN-ID'].values) if val in seq.id]
        if len(num) > 0:
            with open(fasta_rep, 'a') as data_f:
                data_f.write(f'>{seq.id}\n')
                data_f.write(f'{seq.seq}\n')
            with open(fasta_new, 'a') as data_f:
                header = f'{sel_dat["pop"].values[num[0][0]]}|{num[0][1]}'
                data_f.write(f'>{header}\n')
                data_f.write(f'{seq.seq}\n')


    

