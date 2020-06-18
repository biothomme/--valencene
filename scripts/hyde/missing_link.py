#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

this script enables to link SNP information from excel sheets with 
output data from HyDe (after processing with dioloncello.py)


"""
__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'

import os
import pandas as pd
import argparse

# wordir = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
# os.chdir(wordir)
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="excel_file", required=True,
                    help="SNP data with information as exported excel-sheet")
parser.add_argument("-c", action="store", dest="csv_file", required=True,
                    help="output of dioloncello.py as a csv_sheet")

args = parser.parse_args()
excel_file = args.excel_file
csv_file = args.csv_file
# excel_file = '2020_04_vitis_USA.txt'
# csv_file = 'hyde/plot/bonf_.05/2020_04_vitis_USA_hyde_hybr.csv'

indici = []
datae = []
with open(excel_file) as coffee:
    soup = coffee.readlines()   # soup is the list of lines of the file
    for num, spice in enumerate(soup):
        if num <= 4:
            indx = spice.split('\n')[0].split('\t')[0]
            indici.append(indx)
            datx = spice.split('\n')[0].split('\t')[1:]
            datae.append(datx)

reference = pd.DataFrame(data = datae, index = indici).transpose()
reference['cut_names'] = [''.join(eidi.split('_'))[-10:].lower() for eidi in reference['ADN-ID']]
reference.to_csv(f'{excel_file.split(".txt")[0]}_summary.csv')
hybrid_list = pd.read_csv(csv_file)

ind_trip_file = f'{csv_file.split(".csv")[0]}_ind_trp.csv'
with open(ind_trip_file, 'a') as it_f:
    it_f.write(f'parent1,hybrid,parent2,hyb_ind_name,hyb_adn_id,hyb_taxon,hyb_region\n')
for p1, p2, hyb, inds in zip(hybrid_list['p1'], 
                       hybrid_list['p2'], 
                       hybrid_list['hyp_hybrid'],
                       hybrid_list['individuals']):
    inz = inds.replace("[\'", "").replace("\']", "").split("\', \'")
    for ind in inz:
        hybool = reference['cut_names'] == ind
        hynd = reference.loc[hybool, 'Name-long'].values[0]
        hadnid = reference.loc[hybool, 'ADN-ID'].values[0]
        hytax = reference.loc[hybool, 'taxon'].values[0]
        hyloc = reference.loc[hybool, 'region-geo'].values[0]
        with open(ind_trip_file, 'a') as it_f:
            it_f.write(f'{p1},{hyb},{p2},{hynd},{hadnid},{hytax},{hyloc}\n')

triple_list = pd.read_csv(ind_trip_file)

sellist = []
for focal_hyb in set(triple_list['hybrid']):
    order = [x.lower().replace('.', '').replace(' ', '') for x in reference['Pop-celine'].values]
    focal_sel = [i for i, ordo in enumerate(order) if ordo == focal_hyb]
    focal_frame = reference.iloc[focal_sel,:]
    focal_trip = triple_list[triple_list['hybrid'] == focal_hyb]
    for cou, line in focal_frame.iterrows():
        if line['ADN-ID'] not in focal_trip['hyb_adn_id'].values:
            sellist.append(cou)
outliers = reference.iloc[sellist,:]
outliers.to_csv(f'{csv_file.split(".csv")[0]}_not_sign_in_trp.csv')

#print(f'~~~ Finished: {vcf_file} was successcully converted to fasta format: {vcf_file.split(".vcf")[0]}.fasta. ~~~')


