#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  16 13:13:51 2020

this script is made to convert a fastafile to a nexus file to be
used in phylonet (MLE_BiMarkers).
the script produces 01 encoded sequences in nexus format, by 
filtering for biallelic markers and selecting 1 ind per taxon
to be listed afterwards.


"""
__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'

import os
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
import random
from itertools import chain

# makes a subsetted list of individuals within taxa
def indiana_jones(alignment,
                  spset_noasia,
                  spnum):
    indi = []
    stop = 1
    road_kill = 'species_in_set'
    for count, seq in enumerate(alignment):
        eidi = seq.id.split('\\n')[0].split('|')
        eindi = ''.join(eidi[1].split('_'))[-10:].lower()
        eispi = ''.join(eidi[0].split('.')).lower()
        if eispi in spset_noasia:
            if eispi != road_kill:
                road_kill = eispi
                sp_sample = [eindi]
                stop = 1
            else:
                sp_sample.append(eindi)
                stop += 1
            if stop == spnum: 
                indi.append(sp_sample)
        else:
            indi.append([eindi])
    return indi

# produces the final nexus file
def nexico(data_frepname,
       prenex,
       postnex,
       rep_alignment, 
       individual, 
       name_len):
    with open(data_frepname, 'a') as data_f:
        data_f.write(prenex)
    for count, seq in enumerate(rep_alignment):
        eindi = individual[count]
        
        rep_seq = [x for i, x in enumerate(seq.seq) if i not in triallelic]
        bin_seq = ''.join(['1' if rp == st else '0' for rp, st in zip(rep_seq, rep_standard)])
        # eispi = ''.join(eidi[0].split('.')).lower()
        # rewrite the fastafile for PhyloNet as nexus
        with open(data_frepname, 'a') as data_f:
            data_f.write(f'{eindi}{" "*int(name_len + 1 - len(eindi))}')
            data_f.write(f'{bin_seq}\n')
        ## reference all ind to specific populations/taxa
        # with open(map_frepname, 'a') as map_f:
        #     if count+1 != len(alignment):
        #         map_f.write(f'{eindi}\t{eispi}\n')
        #     else:
        #         map_f.write(f'{eindi}\t{eispi}')
    with open(data_frepname, 'a') as data_f:
        data_f.write(postnex)

AMUR = ['vbg176', 'b00f6ps', 'vbg88', 'vbg91']
FLEX = ['vbg184', 'b00erb3']
PIAROM = ['b00er3n', 'vbg267', 'b00er4t']

# os.chdir('/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis')

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="file_name", required=True,
                    help="Name of the input fasta file")

args = parser.parse_args()
file_name = args.file_name

# file_name = '2020_04_vitis_imputed.txt.fasta'

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
    ind.append(eindi)

spnum = 2  # maximal number of ind per species
max_ret = 4

spset_noasia = [grape for grape in set(spac) if grape != 'vitisasia']

for ret in range(max_ret):    
    repsdir = f'{"/".join(file_name.split("/")[:-1])}/{file_name.split("/")[-1].split(".")[0]}_pn_rep_mret{ret}'
    if not os.path.exists(repsdir):
        os.makedirs(repsdir)
    
    print(f'~~~ fasta-file {file_name} was converted to {number_of_replicates} \
    replicates and set for max {ret} reticulation in: {repsdir} ~~~\n ')
        
    # This loop produces all the replicates
    for repo in range(number_of_replicates):
        cou = str(repo + 1)
        lendif = len(str(number_of_replicates)) - len(cou)
        if lendif > 0:
            cou = f'{"0"*lendif}{cou}'
        repdir = f'{repsdir}/rep_{cou}'
        os.makedirs(repdir)
        print(f'  - {repdir}/pn_data.nexus -')
    
        # coordinates of the biallelic loci of this rep within original file (count starts with 0)
        loci_frepname = f'{repdir}/rep_biallelic_loci.txt'
        # nexus_file for phylonet-input
        data_frepname = f'{repdir}/pn_data.nexus'
        ## map file referecing taxon and ind.
        # map_frepname = f'{repdir}/pnw_imap.txt'
        
        # make a replicate
        replicate = list(chain.from_iterable([random.sample(\
            [i for i, x in enumerate(spac) if x == sp], \
            k=int(spnum)) for sp in spset_noasia]))
        # Vitis amurensis (most basal sp in Vitis asia - monophyletic, also cultivated)
        amur = list((random.sample(\
            [i for i, x in enumerate(ind) if x in AMUR], \
            k=1)))
        # Vitis flexuosa (Vitis asia with very different distribution and no big importance as cultivar)
        flex = list((random.sample(\
            [i for i, x in enumerate(ind) if x in FLEX], \
            k=1)))
        # Vitis piasezkii and romanetii (Vitis asia species with reticulate evolution, close to amurensis)
        piarom = list((random.sample(\
            [i for i, x in enumerate(ind) if x in PIAROM], \
            k=1)))
        replicate += amur + flex + piarom
        rep_alignment = MultipleSeqAlignment([alignment[repl] for repl in replicate])
        
        taxa = [x for x in spset_noasia] + ['vamur', 'vflex', 'vpiarom'] # extended taxonset
        
        # which loci are mono- or triallelic?
        triallelic = []
        for base_count in range(len(rep_alignment[0].seq)):
            occ_bases = set(list(rep_alignment[:,base_count]))
            if len([x for x in occ_bases if x.upper() not in ['A','T','G','C']]) > 0:
            	print('You have bases which are not ATGC. This script does not account for N or -.\
please chage the script and rerun!')
            else:
                if len(occ_bases) != 2:
                    triallelic.append(base_count)
        
        print(f'    -> there were {len(triallelic)} mono- or triallelic loci reported. those will be \
removed from the alignment for this replicate. A list of the used loci is stored in: {loci_frepname}')
        loci_kept = [i for i, _ in enumerate(rep_alignment[0].seq) if i not in triallelic]
        with open(loci_frepname, 'a') as loci_f:
            loci_f.write(','.join([str(l) for l in loci_kept]))
        
        # first sequence in alignment is by definition standard and will have only 1 in sequence (reference allele)
        rep_standard = [x for i, x in enumerate(rep_alignment[0].seq) if i not in triallelic]
    
        # obtain length of longest ind name
        indi = indiana_jones(rep_alignment, spset_noasia, spnum)
        taxind_postnex = [','.join(x) for x in indi]
        single_ind = [i for subi in indi for i in subi]
        name_len = max([len(y) for y in single_ind])
        
        # defining nexus preambel and phylonet block
        prenex = f'''#NEXUS\nBegin data;\n\
Dimensions ntax={rep_alignment.__len__()} nchar={len(loci_kept)};\n\
Format datatype=dna symbols="012" missing=? gap=-;\n\
Matrix\n\n'''
        postnex = f''';End;\n\n\
BEGIN PHYLONET;\n\
Nexus_Out \"{repdir}/after_pn.nex\";
MLE_BiMarkers -taxa ({",".join(single_ind)}) -pseudo -mnr 10 -mr {ret} -tm \
<{"; ".join([f"{t}:{i}" for t, i in zip(taxa, taxind_postnex)])}>;\n\
END;'''
    
        # writing the nexus file
        nexico(data_frepname, prenex, postnex, rep_alignment, single_ind, name_len)

print(f'~~~ Done! ~~~\n ')



