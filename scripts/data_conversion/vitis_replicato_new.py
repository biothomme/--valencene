#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 09:00:59 2020

This script is thought to produce replicates of SNP datasets (here: vitis)
by using structure analysis probabilities as well as ibs clusters to avoid
a close realtionships

ATTENTION: there is a bug in the loop! although it works very well, sometimes
replicates are produced, which are not possible to finish. That is because
west and east taxa compete for some shared clusters. this in estimated 1 out 
of 20 cases leads to the run being stuck. apart from that it runs well.

we sample weighted out of the clusters! the weight is the probability, to
which an individual was assigned to a taxon.

@author: Thomsn
"""
# choices for sativa
# dataset: vveast = A3 >= .75
# dataset: vvwest = A3 < .75

__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'

import os
import numpy as np
import pandas as pd
import random
import argparse

TOTAL_REP = 20
NUM_EUR = 5 # how many samples each of vveast, vvwest, vsylveast, vslylvwest
NUM_ASIA = 1
FRESHOLD = .2
VITISASIA = ['V. piasezkii',
             'V. armata',
             'V. davidii',
             'V. flexuosa',
             'V. romanetii',
             'V. balansaeana']

VITISONEASIA = ['V. yeshanensis',
                'V. coignetiae',
                'V. thunbergii',
                'V. pentagona']

VITISUSA = ['V. muscadinia']




parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-y", action="store", dest="infile_sylv", required=True,
                    help="structure excel sheet for sylvestris")
parser.add_argument("-a", action="store", dest="infile_sat", required=True,
                    help="structure excel sheet for sativa")
parser.add_argument("-c", action="store", dest="clusters", required=True,
                    help="file containing all assigned clsuters as csv")
parser.add_argument("-i", action="store", dest="ind_summary", required=True,
                    help="summary csv file with all individuals used")

args = parser.parse_args()
infile_sylv = args.infile_sylv
infile_sat = args.infile_sat
clusters = args.clusters
ind_summary = args.ind_summary

# wordir = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/additional_information/structure'
# infile_sylv = 'Structure-SYLVESTRIS-v1.csv'
# infile_sat = 'Structure-SATIVA-v1.csv'
# clusters = '../ibs_cluster/clusters_vitis_13sativa_10sylvestris.csv'
# ind_summary = '../../2020_04_vitis_USA_summary.csv'

# os.chdir(wordir)

ind_dat = pd.read_csv(ind_summary)
clus_dat = pd.read_csv(clusters)
st1 = ['_'.join(grob.split('GT')[0].split('x')) for grob in clus_dat['adn_id']] # reconvert the x_names
clus_dat['adn_id'] = [f'{grob}GT' if 'Guruchay_' in grob else grob for grob in st1]

rep_path = f'{ind_summary.split(".csv")[0]}_{TOTAL_REP}_reps'
if not os.path.exists(rep_path):
    os.makedirs(rep_path)
## sampling europe ## get replicates for the species v. sylvestris and v. sativa
for rep_num in range(TOTAL_REP):
    repl_list = []
    tclust_list = []
    if len(str(rep_num + 1)) == 1:
        nr = f'0{rep_num + 1}'
    else:
        nr = f'{rep_num + 1}'
    new_path = f'{rep_path}/rep_{nr}'
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    if not os.path.exists(f'{new_path}/replicate_data.csv'):
        total_clusters = []
        ind_sampled = []
        tax_sampled = []
        for spc, species in enumerate([infile_sylv, infile_sat]):
            ## first assign any individual to a group 
            struct_dat = pd.read_csv(species)
            if spc == 0:
                TAXA = ['vsylveast', 'vsylvwest']
            else:
                TAXA = ['vvwest', 'vveast']
            
            after_samp = struct_dat.copy()[['ADN-ID', 'ADN-ID', 'ADN-ID', 'ADN-ID', 'ADN-ID']]
            after_samp.columns = ['ADN-ID', 'taxon', 'pop', 'cluster', 'prob_assigned']
            
            clust_occurrence = pd.DataFrame(data = 0, index = list(range(1, max(clus_dat['ward.D'])+1)), columns = TAXA)
            for _, indi in struct_dat.iterrows():
                if spc == 1:
                    criteria = pd.Series(data = [indi.loc['Structure-A1 BALK'],
                          indi.loc['Structure-A2 WCEU'],
                          indi.loc['Structure-A3 MFEA']],
                          index = ['balk',
                                   'wceu',
                                   'mfea'])
                    sig_grp = 'mfea'
                else:
                    criteria = pd.Series(data = [indi.loc['Structure A1 EMCA'],
                          indi.loc['Structure A2 WCEU']],
                          index = ['emca',
                                   'wceu'])
                    sig_grp = 'wceu'
                save_drop = criteria[criteria > .15]
                new_criteria = save_drop / sum(save_drop)   
                samplate = np.cumsum(new_criteria)
                rnum = random.random()
                choice_int = int([i for i, value in enumerate(samplate) if rnum <= value][0])
                grp = new_criteria.index[choice_int]
                after_samp.loc[after_samp['ADN-ID'] == indi['ADN-ID'], 'pop'] = grp
                after_samp.loc[after_samp['ADN-ID'] == indi['ADN-ID'], 'prob_assigned'] = new_criteria[choice_int]
                selclu = clus_dat.loc[clus_dat['adn_id'] == indi['ADN-ID'], 'ward.D2'].values[0]
                after_samp.loc[after_samp['ADN-ID'] == indi['ADN-ID'], 'cluster'] = selclu
                if grp == sig_grp:
                    after_samp.loc[after_samp['ADN-ID'] == indi['ADN-ID'], 'taxon'] = TAXA[1]
                    clust_occurrence.loc[selclu, TAXA[1]] += 1
                else:
                    after_samp.loc[after_samp['ADN-ID'] == indi['ADN-ID'],'taxon'] = TAXA[0]
                    clust_occurrence.loc[selclu, TAXA[0]] += 1
    
            
            ## now start to sample
            clust_occurrence = clust_occurrence[(clust_occurrence[TAXA[0]] != 0) | (clust_occurrence[TAXA[1]] != 0)]
            clus_stat = pd.DataFrame(data = 0,
                                     columns = [TAXA[0], TAXA[1], 'total'],
                                     index = ['sampled', 'necessary', 'available', 'overlap', 'repeat'])
            
            twin_stat = clus_stat.copy()
            for taxon in TAXA:
                clsa = twin_stat.at['sampled',taxon]
                clne = NUM_EUR - clsa
                clav = sum(clust_occurrence.loc[:,taxon] != 0)
                clov = sum((clust_occurrence.loc[:,TAXA[0]] != 0) & (clust_occurrence.loc[:,TAXA[1]] != 0))
                if clne - clav > 0:
                    clre = clne - clav
                else:
                    clre = 0
                clus_stat.loc[:,taxon] = [clsa, clne, clav, clov, clre]
            twin_stat = clus_stat.copy()
            clsa = sum(twin_stat.iloc[0,0:2])
            clne = sum(twin_stat.iloc[1,0:2])
            clav = sum(twin_stat.iloc[2,0:2]) - clov
            clov = clov
            clre = 0
            clus_stat.loc[:,'total'] = [clsa, clne, clav, clov, clre]
            clusters_sampled = []
            rpc = 0
            while (len(clusters_sampled) < 10) | (rpc > 1000):
                rpc =+1
                if rpc == 1000:
                    print(f'replicate {nr} was hook in the loop and stopped. retry...')
                clust_occ = clust_occurrence.copy().loc[[True if co not in clusters_sampled else False for co in clust_occurrence.index],:]
                the_clus = int(random.sample(list(clust_occ.index), 1)[0])
                # clusters_sampled += [the_clus]
                clust_occ = clust_occ.copy().loc[[True if co != the_clus else False for co in clust_occ.index],:]
                this_sample = after_samp.copy().loc[[True if co == the_clus else False for co in after_samp['cluster'].values],:]
                # HERE!!!
                probass = this_sample.copy()['prob_assigned'] / sum(this_sample.copy()['prob_assigned'])
                samplate = np.cumsum(probass)
                rnum = random.random()
                choice_int = int([i for i, value in enumerate(samplate.values) if rnum <= value][0])
                the_ind = this_sample['ADN-ID'].values[choice_int]
                # the_ind = random.sample(list(this_sample['ADN-ID'].values),1)[0]
                the_tax = this_sample.loc[this_sample['ADN-ID'] == the_ind, 'taxon'].values[0]
                # after_sample = after_samp.copy().loc[[True if co not in clusters_sampled else False for co in after_samp['cluster'].values],:]
                booltot = clus_stat.copy().at['necessary', 'total'] <= clus_stat.copy().at['available', 'total']
                gogo = False
                if the_clus not in clusters_sampled:
                    gogo = True
                elif clus_stat.at['repeat',the_tax] > 0:
                    twin_stat = clus_stat.copy()
                    clre = twin_stat.at['repeat',taxon] - 1
                    lline = list(twin_stat.loc[:,'total'].iloc[0:4].values) + [clre]
                    clus_stat.loc[:,'total'] = lline
                go = False
                if booltot:
                    golist = []
                    for taxon in TAXA:
                        if taxon == the_tax:
                            clsa = clus_stat.copy().at['sampled',taxon] + 1
                        else:
                            clsa = clus_stat.copy().at['sampled',taxon]
                        clne = NUM_EUR - clsa
                        clav = sum(clust_occ.loc[:,taxon] != 0)
                        bool1 = clne <= clav
                        golist.append(bool1)
                    if all(golist):
                        go = True
                if go & gogo:
                    if clus_stat.at['sampled',the_tax] != NUM_EUR:
                        twin_stat = clus_stat.copy()
                        for taxon in TAXA:
                            if taxon == the_tax:
                                clsa = twin_stat.at['sampled',taxon] + 1
                            else:
                                clsa = twin_stat.at['sampled',taxon]
                            clne = NUM_EUR - clsa
                            clav = sum(clust_occ.loc[:,taxon] != 0)
                            clov = sum((clust_occ.loc[:,TAXA[0]] != 0) & (clust_occ.loc[:,TAXA[1]] != 0))
                            clre = twin_stat.at['repeat',taxon]
                            clus_stat.loc[:,taxon] = [clsa, clne, clav, clov, clre]
                        twin_stat = clus_stat.copy()
                        clsa = sum(twin_stat.iloc[0,0:2])
                        clne = sum(twin_stat.iloc[1,0:2])
                        clav = sum(twin_stat.iloc[2,0:2]) - clov
                        clov = clov
                        clre = 0
                        clus_stat.loc[:,'total'] = [clsa, clne, clav, clov, clre]
                        clusters_sampled += [the_clus]
                        ind_sampled += [the_ind]
                        tax_sampled += [the_tax]
            total_clusters += clusters_sampled
    
        ## get replicates for the asian taxa
        ## 1 ind each of V. piasezkii, V. armata, V. davidii, V. flexuosa, V. romanetii, V. balansaeana
        ## 1 ind of V. muscadinia
        ## 1 ind of the group V. yeshanensis, V. coignetiae, V. thunbergii, V. pentagona
        
        ind_dat['taxon'] = [' '.join(tax.split(' ')[0:2]) for tax in ind_dat['taxon']]  #get rid of blank spaces after taxon names
        ind_dat['taxon'] = ['V. balansaeana' if tax == 'V. balanseana' else tax for tax in ind_dat['taxon']] # get rid of two balansaeanas with different spellings
    
    
        for spec in VITISASIA + VITISUSA + [VITISONEASIA]:
            if spec == 'V. muscadinia':
                this_tax = 'vitisusa'
            else:
                this_tax = 'vitisasia'
            if len(spec) == 1:
                sep_dat = ind_dat[ind_dat['taxon'] == spec]
            else:
                sep_bool = [True if sp in spec else False for sp in ind_dat['taxon']]
                sep_dat = ind_dat[sep_bool]
            offer = sep_dat['ADN-ID'].values
            pick = random.sample(list(offer), NUM_ASIA)
            ind_sampled += pick
            total_clusters += ['NA']
            tax_sampled += [this_tax]
        
        ## subset the datasets for the sampled ind
        repl_bool = [True if adn in ind_sampled else False for adn in ind_dat['ADN-ID']]
        repl_dat = ind_dat.copy()[repl_bool]
        clustri = pd.Series(data = total_clusters, index = ind_sampled)
        taxtri = pd.Series(data = tax_sampled, index = ind_sampled)
        repl_dat['cluster'] = clustri[repl_dat['ADN-ID']].values
        repl_dat['pop'] = taxtri[repl_dat['ADN-ID']].values
        repl_dat = repl_dat.drop('Unnamed: 0',1)
        # sorted(ind_sampled)
        repl_dat.to_csv(f'{new_path}/replicate_data.csv')
