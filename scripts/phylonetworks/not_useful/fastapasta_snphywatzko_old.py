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
parser.add_argument("-i", action="store", dest="file_name", required=True,
                    help="Name of the input fasta file")

args = parser.parse_args()
file_name = args.file_name

# file_name = '2020_04_vitis_imputed.txt.fasta'
data_fname = f'{file_name.split(".")[0]}_pnw.phy'
map_fname = f'{file_name.split(".")[0]}__pnw_imap.txt'

number_of_replicates = 20

alignment = AlignIO.read(file_name, 'fasta')
spac = []
ind = []
map_str = ''

for count, seq in enumerate(alignment):
    eidi = seq.id.split('\\n')[0].split('|')
    eindi = ''.join(eidi[1].split('_'))[-10:].lower()
    eispi = ''.join(eidi[0].split('.')).lower()
    spac.append(eispi)


spnum = 10 # maximal number of ind per species
maxnum = len(list(set(spac))) * spnum

repsdir = f'{"/".join(file_name.split("/")[:-1])}/{file_name.split("/")[-1].split(".")[0]}_pnw_rep'
if not os.path.exists(repsdir):
    os.makedirs(repsdir)

print(f'~~~ fasta-file {file_name} was converted to {number_of_replicates} \
replicates in: {repsdir} ~~~\n ')

# This loop produced all the replicates and
for repo in range(number_of_replicates):
    cou = str(repo + 1)
    lendif = len(str(number_of_replicates)) - len(cou)
    if lendif > 0:
        cou = f'{"0"*lendif}{cou}'
    repdir = f'{repsdir}/rep_{cou}'
    os.makedirs(repdir)
    replicate = list(chain.from_iterable([random.sample(\
        [i for i, x in enumerate(spac) if x == sp], \
        k=int(maxnum/len(set(spac)))) for sp in set(spac)]))
    
    spec = []
    ind = []
    
    data_frepname = f'{repdir}/pnw_data.phy'
    map_frepname = f'{repdir}/pnw_imap.txt'
    print(f'  - {repdir}/pnw_data.phy - {repdir}/pnw_imap.txt -')
    with open(data_frepname, 'a') as data_f:
        data_f.write(f'{len(list(set(spac)))} {len(alignment[0].seq)}\n')
        
    with open(map_frepname, 'a') as map_f:
        map_f.write(f'traits\tspecies\n')
    
    for count, seq in enumerate(alignment):
        if count in replicate:
            eidi = seq.id.split('\\n')[0].split('|')
            eindi = ''.join(eidi[1].split('_'))[-10:].lower()
            eispi = ''.join(eidi[0].split('.')).lower()
            spec.append(eispi)
            ind.append(eindi)
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
