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

# number of reticulations
RET = 0
# number of runs
MAX_NUM_RUN = 20

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="rep_dir", required=True,
                    help="Name of the input fasta file")
parser.add_argument('-n', action='store', dest='topo',
                    help='nexus file with topology if given')

args = parser.parse_args()
BIG_DIR = args.rep_dir
TOPO = args.topo

# BIG_DIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonet/2020_04_vitis_USA_summary_20_reps'
os.chdir(BIG_DIR)
# TOPO = True



rep_list = [diro for diro in os.listdir() if 'rep_' in diro]
for repsdir in rep_list:
#    repsdir = f'{rep}'
    file_name1 = f'{repsdir}/rep_new_asgnmt.fasta' # change to std
    file_name2 = f'{repsdir}/rep_new_asgnmt.fasta'
    for file_name in [file_name1, file_name2]:
        alignment = AlignIO.read(file_name, 'fasta')
        spac = []
        ind = []
        map_str = ''
        
        for count, seq in enumerate(alignment):
            eidi = seq.id.split('\\n')[0].split('|')
            eindi = ''.join(eidi[1].split('_'))[-10:].lower()
            eispi = ''.join(''.join(eidi[0].split('.')).split('2')).lower()
            spac.append(eispi)
            ind.append(eindi)
                
        #for ret in range(max_ret + 1):    
        
        
        #    if not os.path.exists(repsdir):
        #        os.makedirs(repsdir)
            
        print(f'~~~ fasta-file {file_name} was set for max {RET} reticulation in: {repsdir} ~~~\n ')
        
        loci_frepname = f'{repsdir}/rep_{file_name.split("_")[-2]}_biallelic_loci.txt'
        # nexus_file for phylonet-input
        data_frepname = f'{repsdir}/rep_{file_name.split("_")[-2]}_data_{RET}_ret.nexus'
        ## map file referecing taxon and ind.
        # map_frepname = f'{repdir}/pnw_imap.txt'
        
        # which loci are mono- or triallelic?
        triallelic = []
        for base_count in range(len(alignment[0].seq)):
            occ_bases = set(list(alignment[:,base_count]))
            if len([x for x in occ_bases if x.upper() not in ['A','T','G','C']]) > 0:
            	print('You have bases which are not ATGC. This script does not account for N or -.\
please change the script and rerun!')
            else:
                if len(occ_bases) != 2:
                    triallelic.append(base_count)
        
        print(f'    -> there were {len(triallelic)} mono- or triallelic loci reported. those will be \
removed from the alignment for this replicate. A list of the used loci is stored in: {loci_frepname}')
        loci_kept = [i for i, _ in enumerate(alignment[0].seq) if i not in triallelic]
        with open(loci_frepname, 'a') as loci_f:
            loci_f.write(','.join([str(l) for l in loci_kept]))
        
        # first sequence in alignment is by definition standard and will have only 1 in sequence (reference allele)
        rep_standard = [x for i, x in enumerate(alignment[0].seq) if i not in triallelic]
        
        # obtain length of longest ind name
        name_len = max([len(y) for y in ind])
        taxa = list(set(spac))
        taxind_postnex = [','.join([ind[i] for i, sp in enumerate(spac) if sp == spe]) for spe in taxa]
        
        
        # defining nexus preambel and phylonet block
        prenex = f'''#NEXUS\nBegin data;\n\
Dimensions ntax={alignment.__len__()} nchar={len(loci_kept)};\n\
Format datatype=dna symbols="012" missing=? gap=-;\n\
Matrix\n\n'''
        ret = 0
        postnex = f''';End;\n\n\
BEGIN PHYLONET;\n\
MLE_BiMarkers -taxa ({",".join(ind)}) -pseudo -mnr 10 -mr {RET} -tm \
<{"; ".join([f"{t}:{i}" for t, i in zip(taxa, taxind_postnex)])}>;\n\
END;'''
        
        # extending postex, if TOPO is given
        if TOPO:
            topo_file = f'{repsdir}/rep_{file_name.split("_")[-2]}_topo.nex'
            with open(topo_file,'r') as tfile:
                soup = tfile.readlines()
                for molecules in soup:
                    if '[&U]' in molecules:
                        atom = molecules.split('[&U] ')[-1]
            postnex = f''';End;\n\n\
BEGIN NETWORKS;\n\
Network net1 = {atom}\
END;\n\n\
BEGIN PHYLONET;\n\
MLE_BiMarkers -taxa ({",".join(ind)}) -pseudo -esptheta -mnr {MAX_NUM_RUN} -mr {RET} -snet (net1) -tm \
<{"; ".join([f"{t}:{i}" for t, i in zip(taxa, taxind_postnex)])}>;\n\
END;'''
        
        else:
            postnex = f''';End;\n\n\
BEGIN PHYLONET;\n\
MLE_BiMarkers -taxa ({",".join(ind)}) -pseudo -mr {RET} -tm \
<{"; ".join([f"{t}:{i}" for t, i in zip(taxa, taxind_postnex)])}>;\n\
END;'''

        # writing the nexus file
        nexico(data_frepname, prenex, postnex, alignment, ind, name_len)
        
print(f'~~~ Done! ~~~\n ')



