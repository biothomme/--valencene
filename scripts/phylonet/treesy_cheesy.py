#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 18:09:40 2020

@author: Thomsn
"""
import argparse
import numpy as np
import os
import pandas as pd


def subtree(tr):
    st =')'.join('('.join(tr.split('(')[1:]).split(')')[:-1])
    return st

def tree_splits(tr):
    comma_split = subtree(tr).split(',')
    opener = np.cumsum([spl.count('(') for spl in comma_split])
    closer = np.cumsum([spl.count(')') for spl in comma_split])
    branching_point = [i for i, (op, cl) in enumerate(zip(opener, closer)) if op <= cl][0] + 1
    return branching_point

def tree_splitter(tr):
    comma_split = subtree(tr).split(',')
    branching_point = tree_splits(tr)
    lefter = ','.join(comma_split[:branching_point])
    righter = ','.join(comma_split[branching_point:])
    return [lefter, righter]

def taximize(subtr):
    TAXA = ['vitisasia',
            'vitisusa',
            'vveast',
            'vvwest',
            'vsylveast',
            'vsylvwest']
    given_tax = [taxon for taxon in TAXA if taxon in subtr]
    return given_tax

def get_the_node(tr):
    import re
    extension = tr.split(')')[-1]
    i_node = re.match(r'I[0-9]*', extension)
    return i_node[0]

def extract_brl(subtr):
    brl = subtr.split(':')[-1]
    return brl

def tree_list(tr):
    jj = tree_splitter(tr)
    trl = [get_the_node(tr),
           taximize(jj[0]),
           extract_brl(jj[0]),
           jj[0],
           taximize(jj[1]),
           extract_brl(jj[1]),
           jj[1],
           tr]
    return trl
    
def tree_sep_branches(tr):
    BRCOL = ['node',
         'branch_1',
         'brl_1',
         'str_1',
         'branch_2',
         'brl_2',
         'str_2',
         'tree_str']
    binaframe = pd.DataFrame(columns = BRCOL)
    
    max_nodes = tr.count(',')
    binaframe = pd.DataFrame(columns = BRCOL)
    binaframe.loc[str(len(binaframe) + 1)] = tree_list(tr)
    while len(binaframe) < max_nodes:
        for i in range(len(binaframe)):
            for side in ['str_1', 'str_2']:
                seltree = binaframe[side].values[i]
                if seltree not in list(binaframe['tree_str'].values):
                    if len(taximize(seltree)) > 1:
                        binaframe.loc[str(len(binaframe) + 1)] = tree_list(seltree)
    return binaframe

def compare_topos(tree_csv, topology=None):
    if topology:
        ts_ref = topology
    else:    
        treeset = tree_csv['topology'].values
        ts_ref = treeset[0]
    reference = tree_sep_branches(ts_ref).loc[:,['node','branch_1','branch_2']]
    for co, trs in tree_csv.iterrows():
        tro = tree_sep_branches(trs['topology'])
        share_node = []
        for i in range(len(reference)):
            curref = reference.copy().iloc[i,:]
            for (rb, lb) in [('branch_1','branch_2'), ('branch_2','branch_1')]:
                res = tro[rb].copy().values
                les = tro[lb].copy().values
                lefters = [j \
                       for j, (lf, rg) in \
                           enumerate(zip(res,les)) \
                               if (lf == curref['branch_1'] and rg == curref['branch_2'])]
                if len(lefters) > 0:
                    share_node.append(tro['node'].copy().values[lefters[0]])
                    break
                elif (rb, lb) == ('branch_2','branch_1'):
                    share_node.append(None)
        reference.loc[:,trs[0]] = share_node
    return reference

def summarize_topos(topology_matrix):
    sel_mat = topology_matrix.iloc[:,3:].copy()
    top_summary = topology_matrix.iloc[:,0:1].copy()
    topo_count = len(sel_mat.iloc[1,:])
    occ_list = []
    for co, toprow in sel_mat.iterrows():
        occs = sum([1 for el in toprow if el != None])
        occ_list.append(occs)
    top_summary.loc[:,'occurences'] = occ_list
    top_summary.loc[:,'percentage'] = np.divide(occ_list, topo_count)
    return top_summary
            
def make_bootstrap_str(tree_csv, topology_matrix = pd.DataFrame()):
    if topology_matrix.empty:
        topology_matrix = summarize_topos(compare_topos(tree_csv))
    template_str = tree_csv['topology'].values[0]
    for co, toporow in topology_matrix.iterrows():
        nod = toporow['node']
        val = f'{toporow["percentage"]:.3f}'
        template_str = val.join(template_str.split(nod))
    return template_str

def check_for_topo(tree_csv, topo_candidate):
    for co, trs in tree_csv.iterrows():
        tro = tree_sep_branches(trs['topology'])
        share_node = []
        for i in range(len(topo_candidate)):
            curref = topo_candidate.copy().iloc[i,:]
            for (rb, lb) in [('branch_1','branch_2'), ('branch_2','branch_1')]:
                res = tro[rb].copy().values
                les = tro[lb].copy().values
                lefters = [j \
                       for j, (lf, rg) in \
                           enumerate(zip(res,les)) \
                               if (lf == curref['branch_1'] and rg == curref['branch_2'])]
                if len(lefters) > 0:
                    share_node.append(tro['node'].copy().values[lefters[0]])
                    break
                elif (rb, lb) == ('branch_2','branch_1'):
                    share_node.append(None)
        topo_candidate.loc[:,trs[0]] = share_node
    return topo_candidate
        
        
def get_all_topos(tree_csv):
    tree_csv = all_trees
    treeset = tree_csv['topology'].values
    total_topoframe = pd.DataFrame()
    if total_topoframe.size < 1:
        nextr = 0
        rank = [1]
        entrance = True
    while entrance:
        reference = tree_sep_branches(treeset[nextr]).loc[:,['node','branch_1','branch_2']]
        reference['topo_rank'] = reference.shape[0] * rank
        reference['best_topo'] = reference.shape[0] * [treeset[nextr].split('\n')[0]]
        if total_topoframe.size < 1:
            total_topoframe = check_for_topo(tree_csv, reference)
        else:
            total_topoframe = pd.concat([total_topoframe,check_for_topo(tree_csv, reference)])
        all_ranks = set(total_topoframe['topo_rank'])
        bool_frame = pd.DataFrame(columns=all_ranks)
        for rk in all_ranks:
            frame_sel = total_topoframe.loc[total_topoframe['topo_rank'] == rk,:].iloc[:,5:]
            bool_frame[rk] = [True if None in frame_sel[coln].values else False for i, coln in enumerate(frame_sel.columns)]
        disc_topologies = [i for i in bool_frame.index if all(bool_frame.iloc[i,:].values) == True]
        rank = [rank[0] + 1]
        if len(disc_topologies) > 0:
            nextr = disc_topologies[0]
        else:
            entrance = False
    return total_topoframe

def all_topos_nodlab(tree_csv, total_topoframe):
    sel_topo_sum = summarize_topos(total_topoframe).copy()
    sel_topo_sum.columns = ['node', 'occurences', 'percentage']
    sbt = make_bootstrap_str(tree_csv, sel_topo_sum)
    return sbt

def inverse_bool(boolean):
    return not boolean

def obt_likelihoods(tree_csv, total_topoframe):
    NEWCOLS = ['rank', 'pseudo_likelihoods', 'replicates', 'best_topology']
    all_ranks = set(total_topoframe['topo_rank'])
    sum_frame = pd.DataFrame(columns = NEWCOLS)
    for rk in all_ranks:
        frame_full = total_topoframe.loc[total_topoframe['topo_rank'] == rk,:]
        frame_sel = frame_full.copy().iloc[:,5:]
        bf = [False if None in frame_sel[coln].values else True for i, coln in enumerate(frame_sel.columns)]
        rf = [coln for coln in frame_sel.columns if not None in frame_sel[coln].values]
        frame_sep = pd.concat([frame_full.copy().iloc[:,0:5], 
                               frame_sel.loc[:, bf], 
                               frame_sel.loc[:, list(map(inverse_bool, bf))]],
                              axis=1)
        bto = all_topos_nodlab(tree_csv.loc[bf,:], frame_sep)
        lk = list(tree_csv.loc[bf,'pseudo_likelihood'].values)
        sum_frame.loc[rk, :] = [rk, lk, rf, bto]
    return sum_frame

def dom_2_topo(tree_csv):
    import re
    TOP_TOPO = '((((vveast,vvwest)I4,vsylveast)I3,vsylvwest)I2,(vitisasia,vitisusa)I1)I0;\n'
    BRANCH_LENGTH_PART = .25
    reference = compare_topos(tree_csv, topology=TOP_TOPO).iloc[:,3:]
    bf = [False if None in reference[coln].values else True for i, coln in enumerate(reference.columns)]
    if any(bf):
        ip_topo = tree_csv.loc[bf,:].iloc[0,:]['topology']
        vsylv = re.search(r'vsylvwest:[0-9]*\.[0-9]*E?\-?[0-9]*',
                          ip_topo).group(0)
        vv = re.search(r'vvwest:[0-9]*\.[0-9]*E?\-?[0-9]*',
                          ip_topo).group(0)
        vvbl = float(re.search(r'[0-9]*\.[0-9]*E?\-?[0-9]*', vv).group(0))
        newvvbls = [i * vvbl for i in [BRANCH_LENGTH_PART, 1-BRANCH_LENGTH_PART]]
        vsylvbl = float(re.search(r'[0-9]*\.[0-9]*E?\-?[0-9]*', vsylv).group(0)) - newvvbls[1]
        sylvec = re.sub(vsylv, f'(vsylvwest:{newvvbls[1]})#H1:{vsylvbl}', ip_topo)
        input_topo = re.sub(vv, f'(vvwest:{newvvbls[1]})#H1:{newvvbls[0]}', sylvec)
    else:
        input_topo = tree_csv.loc[0,:]['topology']
    return input_topo


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="bigdir", required=True,
                    help="big dir with repdirs")
args = parser.parse_args()
BIG_DIR = args.bigdir
NEWCOLS = ['rank', 'pseudo_likelihoods', 'replicates', 'best_topology']

#BIG_DIR = 'Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonet/2020_04_vitis_USA_2_summary_20_reps/'
os.chdir(BIG_DIR)

rep_list = [diro for diro in os.listdir() if 'rep_' in diro]
all_stars = pd.DataFrame(columns = NEWCOLS)
for repsdir in rep_list:
    all_trees = pd.read_csv(f'{repsdir}/tree_sum.csv')
    # topo_mat = compare_topos(all_trees)
    # topo_sum = summarize_topos(topo_mat)
    # bt = make_bootstrap_str(all_trees, topo_sum)
    # best_trees = f'{repsdir}/best_tree_rep.new'
    # with open(best_trees, 'w') as bff:
    #     bff.write(bt)
    tree_csv = all_trees.copy()
    summary_frame = obt_likelihoods(all_trees, get_all_topos(all_trees))
    summary_frame.to_csv(f'{repsdir}/net_0_summary_seedrep.csv')
    all_stars.loc[repsdir,:] = summary_frame.iloc[0,:]
    topo_dir = f'{repsdir}/topologies_new_ret1'
    if not os.path.exists(topo_dir):
        os.makedirs(topo_dir)
    input_topofile = f'{topo_dir}/ret_01_topo.new'
    with open(input_topofile, 'w') as itf:
        itf.write(dom_2_topo(tree_csv))
    
all_stars.to_csv('net_0_bt_samrep.csv')
