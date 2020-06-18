#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  14 13:13:51 2020

this script is made to convert a fastafile to a phylip file to be
used in SNPs2CF and so on phylonetworks.
the script produces replicates with 10 ind per species per rep.
this is neccessary for SNPs2CF


"""
__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'

import os
import argparse
from Bio import AlignIO
import pandas as pd
import random
from itertools import chain


# os.chdir('/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis')

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="bigdir", required=True,
                    help="name of bigdir")

args = parser.parse_args()
BIGDIR = args.bigdir
TYPE_AN = ['std'] # ['new', 'std']
os.chdir(BIGDIR)

rep_list = [diro for diro in os.listdir() if 'rep_' in diro]
for repsdir in rep_list:
    for antype in TYPE_AN:
        file_name = f'{repsdir}/rep_{antype}_asgnmt.fasta'
        # file_name = '2020_04_vitis_imputed.txt.fasta'
        data_fname = f'{file_name.split(".")[0]}_pnw.phy' # {antype}
        map_fname = f'{file_name.split(".")[0]}__pnw_imap.txt' # {antype}
        
        
        alignment = AlignIO.read(file_name, 'fasta')
        spac = []
        ind = []
        map_str = ''
        
        for count, seq in enumerate(alignment):
            eidi = seq.id.split('\\n')[0].split('|')
            eindi = ''.join(eidi[1].split('_'))[-10:].lower()
            eispi = ''.join(eidi[0].split('.')).lower().split('2')[0]
            spac.append(eispi)
            ind.append(eindi)

        print(f'~~~ fasta-file {file_name} was converted ~~~\n ')
        
    
        data_frepname = f'{repsdir}/pnw_data.phy'  # {antype}
        map_frepname = f'{repsdir}/pnw_imap.txt' # {antype}
        print(f'  - {repsdir}/pnw_data.phy - {repsdir}/pnw_imap.txt -')
        with open(data_frepname, 'a') as data_f:
            data_f.write(f'{len(list(set(ind)))} {len(alignment[0].seq)}\n')
            
        with open(map_frepname, 'a') as map_f:
            map_f.write(f'traits\tspecies\n')
        
        for count, seq in enumerate(alignment):
            eidi = seq.id.split('\\n')[0].split('|')
            eindi = ''.join(eidi[1].split('_'))[-10:].lower()
            eispi = ''.join(eidi[0].split('.')).lower()
    
            # rewrite the fastafile for HyDe as phylip - names need to be shortened to last 10 chars
            with open(data_frepname, 'a') as data_f:
                data_f.write(f'{eindi}{" "*int(11 - len(eindi))}')
                data_f.write(f'{seq.seq}\n')
            # reference all ind to specific populations/taxa
            with open(map_frepname, 'a') as map_f:
                if count+1 != len(alignment):
                    map_f.write(f'{eindi}\t{eispi}\n')
                else:
                    map_f.write(f'{eindi}\t{eispi}')
        print(f'~~~ Done! ~~~\n ')
