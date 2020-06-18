#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 13:13:51 2020

this script is made to convert a fastafile to feasible input files (data, map and triples)
for an individual_hyde.py analysis


"""
__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'

# import os
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
data_fname = f'{file_name.split(".")[0]}-data.txt'
map_fname = f'{file_name.split(".")[0]}-map.txt'
triples_fname = f'{file_name.split(".")[0]}-triples.txt'

alignment = AlignIO.read(file_name, 'fasta')
spec = []
ind = []
map_str = ''
for count, seq in enumerate(alignment):
    eidi = seq.id.split('\\n')[0].split('|')
    eindi = ''.join(eidi[1].split('_'))[-10:].lower()
    eispi = ''.join(eidi[0].split('.')).lower()
    spec.append(eispi)
    ind.append(eindi)
    # rewrite the fastafile for HyDe as phylip - names need to be shortened to last 10 chars
    with open(data_fname, 'a') as data_f:

        data_f.write(f'{eindi}{" "*int(11 - len(eindi))}')
        data_f.write(f'{seq.seq}\n')
    # reference all ind to specific populations/taxa
    with open(map_fname, 'a') as map_f:
        if count+1 != len(alignment):
            map_f.write(f'{eindi}\t{eispi}\n')
        else:
            map_f.write(f'{eindi}\t{eispi}')

if len(set(ind)) != len([species[-10:] for species in list(set(ind))]):
    print('Your species ids are too long, please change it and do not trust the produced files!!!')
[f'{x}: {len([i for i in spec if i == x])}' for x in set(spec)]
# finally permutate all taxa to receive all triples. remove outgroup!
outgroup = 'vitisasia'
spelist = [tax for tax in set(spec) if tax != outgroup]
outlist = []
tri_str = ''
for spo in spelist:
    spolist = [sp for sp in spelist if sp != spo]
    for spu in spolist:
        spulist = [sp for sp in spolist if sp != spu]
        spulist
        for spa in spulist:
            if spa not in outlist:
                tri_str += f'{spo}\t{spu}\t{spa}\n'
    outlist.append(spo)
    
with open(triples_fname, 'a') as tri_f:
    tri_f.write(tri_str[:-1])
    

# HyDe showed problems with too many individuals. Therefore smaller replicates
# can be produced: 
# maxnum = 25 # maximal number of ind
# replicate = list(chain.from_iterable([random.sample([i for i, x in enumerate(spec) if x == sp], k=int(maxnum/len(set(spec)))) for sp in set(spec)]))

# data_frepname = f'{file_name.split(".")[0]}-data-rep.txt'
# map_frepname = f'{file_name.split(".")[0]}-map-rep.txt'

# for count, seq in enumerate(alignment):
#     if count in replicate:
#         eidi = seq.id.split('\\n')[0].split('|')
#         eindi = ''.join(eidi[1].split('_'))[-10:].lower()
#         eispi = ''.join(eidi[0].split('.')).lower()
#         spec.append(eispi)
#         ind.append(eindi)
#         # rewrite the fastafile for HyDe as phylip - names need to be shortened to last 10 chars
#         with open(data_frepname, 'a') as data_f:
#             data_f.write(f'{eindi}{" "*int(11 - len(eindi))}')
#             data_f.write(f'{seq.seq}\n')
#         # reference all ind to specific populations/taxa
#         with open(map_frepname, 'a') as map_f:
#             if count+1 != len(alignment):
#                 map_f.write(f'{eindi}\t{eispi}\n')
#             else:
#                 map_f.write(f'{eindi}\t{eispi}')

    
print(f'~~~ fasta-file {file_name} was converted to: ~~~\n \
- {data_fname}\n - {map_fname}\n - {triples_fname}') # \n - {data_frepname}\n - {map_frepname}')
