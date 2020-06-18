#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

this script enables to obtain a distance matrix for pairwise IBS data
in the stucture of an excel sheet (as csv). additionally it saves the
location information of the individuals, which was given in the excel
sheet.


"""
__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'


import os
import pandas as pd
import numpy as np
from itertools import combinations
import argparse

# WORDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/big_data/vitis'
# IBS = 'IBS0-IBS-dist-SATIVA.csv'
# IBS = 'IBS0-IBS-dist-SYLVESTRIS.csv'
# os.chdir(WORDIR)

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="ibs_file", required=True,
                    help="IBS data with information as excel-sheet (csv)")

args = parser.parse_args()
IBS = args.ibs_file

ibs_outf = f'{"/".join(IBS.split("/")[0:-1])}/ibs_mat_{IBS.split("-")[-1].split(".")[0].lower()}.csv'
loc_outf = f'{"/".join(IBS.split("/")[0:-1])}/loc_mat_{IBS.split("-")[-1].split(".")[0].lower()}.csv'

ib_dat = pd.read_csv(IBS)
all_ind = ib_dat.loc[:, 'IID1'].values
all_ind = np.append(all_ind, ib_dat.loc[:, 'IID2'].values)
all_individuals = list(set(all_ind))

ibs_mat = pd.DataFrame(0, index=all_individuals, columns=all_individuals)

for i1, i2 in combinations(all_individuals, 2):
    i1_bool1 = ib_dat.loc[:, 'IID1'] == i1
    bool1_dat = ib_dat.loc[i1_bool1, :]
    i2_bool = bool1_dat.loc[:, 'IID2'] == i2
    if sum(i2_bool) == 0:
        i1_bool2 = ib_dat.loc[:, 'IID2'] == i1
        bool1_dat = ib_dat.loc[i1_bool2, :]
        i2_bool = ib_dat.loc[:, 'IID1'] == i2
    if sum(i2_bool) != 0:
        bool11_dat = bool1_dat.loc[i2_bool, :]
        ibs_values = bool11_dat.loc[:, 'IBS'].values
        ibs_mat.loc[i1, i2] = float(ibs_values)
        ibs_mat.loc[i2, i1] = float(ibs_values)
    


ibs_new = 1 - ibs_mat
loc_dat = pd.DataFrame(columns = ['adn_id', 'country', 'name'])

for i1 in all_individuals:
    ibs_new.loc[i1, i1] = 0
    ind_bool = ib_dat.loc[:, 'IID1'] == i1
    meta = ib_dat.loc[ind_bool, :]
    if len(meta) > 1:
        meta = meta.iloc[0, :]
        if 'uruc' in i1:
            country = 'AZE'  # the guruchay individuals still have no country - google says AZE...
            name = i1
        else:
            country = meta['Country1']
            name = meta['Name1']
    else:
        ind_bool = ib_dat.loc[:, 'IID2'] == i1
        meta = ib_dat.loc[ind_bool, :].iloc[0, :]
        if 'uruc' in i1:
            country = 'AZE'  # the guruchay individuals still have no country - google says AZE...
            name = i1
        else:
            country = meta['Country2']
            name = meta['Name2']
    loc_dat.loc[i1,:] = [i1, country, name]

loc_dat.index = range(len(loc_dat))

ibs_new.to_csv(ibs_outf)
loc_dat.to_csv(loc_outf)
        