#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 09:59:35 2020

this file should visualize differences between the replicates

@author: Thomsn
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re

def standardsort(bigdir, xaxis, sortlist=None, pnw_topo=None):
    os.chdir(bigdir)
    rep_list = [diro for diro in os.listdir() if 'rep_' in diro]
    rep_list.sort()
    if sortlist:
        rep_list = [co for _, co in sorted(zip(sortlist[:-1], rep_list))]        
    rep_mat = pd.DataFrame(columns = xaxis)# np.append(xaxis, 'topo'))
    rep_arr = np.array([])
    
    topolist = []
    tclist = []
    shortopolist = []
    if 'pnw_topo' in locals():
        pntopolist = []
        pntclist = []
        pnshortopolist = []
    for repsdir in rep_list:
        rep_ind = pd.read_csv(f'{repsdir}/replicate_data.csv')
        
        rep_topo = pd.read_csv(f'{repsdir}/net_0_summary_seedrep.csv')
        rtop = rep_topo['best_topology'].iat[0]
        if 'pnw_topo' in locals():
            rep_num = int(repsdir.split('rep_')[1])
            ptop = pnw_topo[rep_num-1]
            apstop = rearrange_topo(ptop)
            pntopolist.append(apstop)
            pntclist.append(TOPOCOL[apstop])
            pnshortopolist.append(TOPOSHORT[apstop])
        abstop = rearrange_topo(rtop)
        topolist.append(abstop)
        tclist.append(TOPOCOL[abstop])
        shortopolist.append(TOPOSHORT[abstop])
        ind_ass = []
        for ind in xaxis:
            if ind in list(rep_ind['ADN-ID'].values):
                indpop = rep_ind.loc[rep_ind['ADN-ID'] == ind, 'pop'].values[0]
                ind_ass.append(float(COL_DIC[indpop]))
            else:
                ind_ass.append(float(0))
        # ind_ass.append(float(TOPOCOL[abstop]))
        
        rep_mat.loc[repsdir,:] = ind_ass
        if len(rep_arr) == 0:
            rep_arr = np.array(ind_ass)
        else:
            rep_arr = np.vstack((rep_arr, np.array(ind_ass)))
    reflist = [float(COL_DIC[''.join(''.join(pop.lower().split(' ')).split('.'))])\
               for pop in all_ind['Pop-celine']]
    rep_mat.loc['reference',:] = reflist
    rep_arr = np.vstack((rep_arr, np.array(reflist)))
    tclist.append(0)
    shortopolist.append('phylonet')
    if 'pnw_topo' in locals():
        pntclist.append(0)
        pnshortopolist.append('phylonetworks')
        plt_rep(rep_mat, rep_arr, rep_list, tclist, shortopolist, pntclist, pnshortopolist)
    else:
        plt_rep(rep_mat, rep_arr, rep_list, tclist, shortopolist)
    return pntclist
    

## change colormap ##
def make_new_cm(com, percent_black):
    from matplotlib import cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    viridis = cm.get_cmap(com, 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    wht = np.array([0/256, 0/256, 0/256, percent_black])
    newcolors[:5, :] = wht
    return ListedColormap(newcolors)


def plt_rep(rep_mat, 
            rep_arr, 
            rep_list, 
            topocollist, 
            toposhortlist, 
            popocollist=None, 
            poposhortlist=None):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    newcmp = make_new_cm('hsv', 0.4)
    newcmp2 = make_new_cm('Set1', 0.0)
    libo = [i for i, collo in enumerate(rep_mat.columns) if sum(rep_mat[collo] != 0) > 1]
    fig = plt.figure(figsize=(40, 6))
    
    ax1 = plt.subplot(111)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="15%", pad=0.15)
    c = ax1.pcolor(rep_arr[:,libo], edgecolors='k', linewidths=.1, cmap = newcmp)
    p = cax.pcolor(np.transpose(np.array([topocollist,topocollist,topocollist])),
                   cmap = newcmp2, edgecolors='k', linewidths=0)
    if (('popocollist' in locals()) & ('poposhortlist' in locals())):
        wax = divider.append_axes("right", size="15%", pad=0.15)
        w = wax.pcolor(np.transpose(np.array([popocollist,popocollist,popocollist])),
                   cmap = newcmp2, edgecolors='k', linewidths=0)
    ax1.set_title('replicates')
    ax1.set_xticks(np.arange(len(libo)) + 0.5)
    ax1.set_yticks(np.arange(len(rep_mat.index)) + 0.5)
    ax1.invert_yaxis()
    cax.invert_yaxis()
    cax.set_xticks([])
    cax.set_yticks([])
    for i in range(len(rep_mat.index)):
        text = cax.text(1.5, i + 0.5, toposhortlist[i],
                       ha="center", va="center", color="k")
        
    if (('popocollist' in locals()) & ('poposhortlist' in locals())):
        wax.invert_yaxis()
        wax.set_xticks([])
        wax.set_yticks([])
        for i in range(len(rep_mat.index)):
            text = wax.text(1.5, i + 0.5, poposhortlist[i],
                           ha="center", va="center", color="k")
    ax1.set_xticklabels(rep_mat.columns[libo], fontsize=5)
    ax1.set_yticklabels(rep_mat.index)

    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    plt.savefig('tree_behaviour_replicates.pdf', bbox_inches='tight')

    
def taximize(subtr):
    TAXA = ['vitisasia',
            'vitisusa',
            'vveast',
            'vvwest',
            'vsylveast',
            'vsylvwest']
    given_tax = [taxon for taxon in TAXA if taxon in subtr]
    return given_tax


def rearrange_topo(given_topo):
    import re
    abs_topo = re.sub(r':?[0-9]*\.[0-9]*E?-?[0-9]*', '', given_topo)
    abs_topo = abs_topo.split('\n')[0]
    openers = abs_topo.count('(')
    for op in range(openers):
        jack = '('.join(abs_topo.split('(')[op+1:])
        commasplit = jack.split(',')
        opi = np.cumsum([el.count('(') for el in commasplit])
        clo = np.cumsum([el.count(')') for el in commasplit])
        judge = [i for i, (op, cl) in enumerate(zip(opi,clo)) if op <= cl][0]
        left = ','.join(commasplit[:judge+1])
        right = ','.join(commasplit[judge+1:])
        rommasplit = right.split(',')
        ropi = np.cumsum([el.count('(') for el in rommasplit])
        rclo = np.cumsum([el.count(')') for el in rommasplit])
        rjudge = [i for i, (op, cl) in enumerate(zip(ropi,rclo)) if op <= cl][0]
        right = ')'.join(right.split(')')[:rjudge+1])
        if left.count(',') < right.count(','):
            abs_topo = 'SPLIIIIT'.join(abs_topo.split(left))
            abs_topo = f'{left}'.join(abs_topo.split(right))
            abs_topo = f'{right}'.join(abs_topo.split('SPLIIIIT'))
        elif left.count(',') == right.count(','):
            let = taximize(left)
            rit = taximize(right)
            taxsort = list(set(let + rit))
            taxsort.sort()
            side = ['right' if tax in rit else 'left' for tax in taxsort][0]
            if side == 'right':
                abs_topo = 'SPLIIIIT'.join(abs_topo.split(left))
                abs_topo = f'{left}'.join(abs_topo.split(right))
                abs_topo = f'{right}'.join(abs_topo.split('SPLIIIIT'))
    return abs_topo

def import_pnwtopo(topofile):
    topo_frame = pd.DataFrame(columns = ['rep', 'topo'])
    repl = None
    with open(topofile, 'r') as tof:
        soup = tof.readlines()
        for molecule in soup:
            if 'rep' in molecule:
                repl = int(molecule.split('rep_')[-1].split('\n')[0])
            else:
                topl = molecule.split('\n')[0]
                topo_frame.loc[len(topo_frame),:] = [repl, topl]
    return topo_frame
            
def commasort(topolist):
    elnum = [num.count(',')+1 for num in sorted(topolist, reverse=True)]
    return [elle for _, elle in sorted(zip(elnum, sorted(topolist, reverse=True)), reverse=True)]

def polytomagic(given_topo):
    import re
    if given_topo.count(',') == given_topo.count(')'):
        print('Topology is already binary.')
        return given_topo
    else:
        abs_topo = re.sub(r':?[0-9]*\.[0-9]*E?-?[0-9]*', '', given_topo)
        abs_topo = abs_topo.split('\n')[0]
        new_topo = f'{abs_topo}'
        openers = abs_topo.count('(')
        elements_series = pd.Series()
        po = -1
        for op in range(openers):
            jack = '('.join(abs_topo.split('(')[op+1:])
            jock = ')'.join(jack.split(')')[:po])
            commers = jack.split(',')
            opi = np.cumsum([el.count('(') for el in commers])
            clo = np.cumsum([el.count(')') for el in commers])
            for o, c in zip(opi, clo):
                if o < c:
                    break
                elif (c > 0 & o == 0):
                    po += 1
            commors = jock.split(',')
            print(jock)
            print(opi)
            print(clo)
            elements = []
            prel = None
            po -=1
            for el in commors:
                if prel:
                    whel = ','.join([prel, el])
                else:
                    whel = el
                if whel.count('(') <= whel.count(')'):
                    elements.append(whel)
                    prel = None
                else:
                    prel = whel
            elements_series.loc[op] = elements
        order = [sum([len([indx for indx in COL_DIC.keys() if indx in elli]) for elli in ellist]) for ellist in elements_series]
        elements_df = pd.DataFrame({'data': elements_series, 'rank': order}).sort_values(by = ['rank'], ascending=False)
        for rel in elements_df['data']:
            old = ','.join(rel)
            new = ','.join(commasort(rel))
            print(old)
            print(new)
            new_topo = new.join(new_topo.split(old))
    return new_topo

def topolsum(topol):
    topol = re.sub(r':?[0-9]*\.[0-9]*E?-?[0-9]*', '', topol)
    levels = topol.split('(')
    topologies = pd.Series()
    for i in range(1, len(levels)):
        recovertree = '('.join(levels[i:])
        commens = recovertree.split(',')
        
        opi = np.cumsum([sp.count('(') for sp in commens])
        clo = np.cumsum([sp.count(')') for sp in commens])
        diffo = [op-cl for op, cl in zip(opi, clo)]
        
        elements = []
        whole = []
        for j, com in enumerate(commens):
            if diffo[j] < 0:
                clom = com.split(')')
                crom = ')'.join(clom[:diffo[j]])
                whole.append(crom)
                wstr = ','.join(whole)
                elements.append(wstr)
                break
            elif diffo[j] == 0:
                whole.append(com)
                wstr = ','.join(whole)
                elements.append(wstr)
                whole = []
            else:
                whole.append(com)
        topologies.loc[i] = elements
    return topologies

def topologic(topol, COL_DIC):
    topologies = topolsum(topol)
    topol = re.sub(r':?[0-9]*\.[0-9]*E?-?[0-9]*', '', topol)
    order = [sum([len([indx for indx in COL_DIC.keys() if indx in elli]) for elli in ellist]) for ellist in topologies]
    elements_df = pd.DataFrame({'data': topologies, 'rank': order}).sort_values(by = ['rank'], ascending=False)
    new_topo = topol
    for rel in elements_df['data']:
        old = ','.join(rel)
        new = ','.join(commasort(rel))
        new_topo = new.join(new_topo.split(old))
    return new_topo

def myparents(mystr, COL_DIC):
    if len([indx for indx in COL_DIC.keys() if indx in mystr]) > 1:
        if not ((mystr[0] == '(') & (mystr[-1] == ')')):
            mystr = f'({mystr})'
        co = 0
        while ((mystr[co] == '(' ) & (mystr[-(co+1)] == ')')):
            co += 1
        if co > 0:
            if mystr.count('),(') >= co:
                mystr = f'({mystr})'
    return mystr

def brootal(topol, COL_DIC):
    outgroup = ['vitisasia','vitisusa']
    topologies = topolsum(topol)
    topol = re.sub(r':?[0-9]*\.[0-9]*E?-?[0-9]*', '', topol)
    total_tax = [indx for indx in COL_DIC.keys() if indx in topol]
    instr = None
    total = None
    searcher = outgroup
    while not total:
        for top in topologies:  # should be sorted by descending length!
            top_tax = [[indx for indx in COL_DIC.keys() if indx in tp] for tp in top]
            top_tax = [tt for tota in top_tax for tt in tota]
            ogs = []
            for i, to in enumerate(top):
                to_tax = [indx for indx in COL_DIC.keys() if indx in to]
                if all([True if og in searcher else False for og in to_tax]):
                    ogs.append(i)
            if ogs:
                if not instr:
                    instr = ','.join([too for e, too in enumerate(top) if e not in ogs])
                    instr = myparents(instr, COL_DIC)
                    outstr = ','.join([too for e, too in enumerate(top) if e in ogs])
                    outstr = myparents(outstr, COL_DIC)
                    if len(top_tax) == len(total_tax):
                        total = f'({outstr},{instr});'
                        break
                    else:
                        searcher = top_tax
                else:
                    if all([True if og in top_tax else False for og in searcher]):
                        iistr = ','.join([too for e, too in enumerate(top) if e not in ogs])
                        iistr = myparents(iistr, COL_DIC)
                        instr = f'{iistr},{instr}'
                        instr = myparents(instr, COL_DIC)

                        if len(top_tax) == len(total_tax):
                            total = f'({outstr},{instr});'
                            break
                        else:
                            searcher = top_tax
    return total

SUMDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
COL_DIC = {'vitisusa':1,
 'vitisusa2':1,
 'vitisasia':2,
 'vveast':3,
 'vvwest':4,
 'vsylveast':5,
 'vsylvwest':6}

TOPOCOL = {'((((vsylveast,vsylvwest),vvwest),vveast),(vitisasia,vitisusa));': 1,
 '((((vsylveast,vveast),vvwest),vsylvwest),(vitisasia,vitisusa));': 2,
 '((((vsylveast,vvwest),vsylvwest),vveast),(vitisasia,vitisusa));': 6,
 '((((vveast,vvwest),vsylveast),vsylvwest),(vitisasia,vitisusa));': 7,
 '(((vsylveast,vveast),(vsylvwest,vvwest)),(vitisasia,vitisusa));': 5,
 '(((vsylveast,vsylvwest),(vveast,vvwest)),(vitisasia,vitisusa));': 3,
 '((((vsylvwest,vvwest),vsylveast),vveast),(vitisasia,vitisusa));': 4}

TOPOSHORT = {'((((vsylveast,vsylvwest),vvwest),vveast),(vitisasia,vitisusa));': '(vas,vu2),(vve,(vvw,(vse,vsw)))',
 '((((vsylveast,vveast),vvwest),vsylvwest),(vitisasia,vitisusa));': '(vas,vu2),(vsw,(vvw,(vse,vve)))',
 '((((vsylveast,vvwest),vsylvwest),vveast),(vitisasia,vitisusa));': '(vas,vu2),(vve,(vsw,(vse,vvw)))',
 '((((vveast,vvwest),vsylveast),vsylvwest),(vitisasia,vitisusa));': '(vas,vu2),(vsw,(vse,(vve,vvw)))',
 '(((vsylveast,vveast),(vsylvwest,vvwest)),(vitisasia,vitisusa));': 'scen_class',
 '(((vsylveast,vsylvwest),(vveast,vvwest)),(vitisasia,vitisusa));': '(vas,vu2),((vse,vsw),(vve,vvw))',
 '((((vsylvwest,vvwest),vsylveast),vveast),(vitisasia,vitisusa));': '(vas,vu2),(vve,(vse,(vsw,vvw)))' }


os.chdir(SUMDIR)
all_ind = pd.read_csv('2020_04_vitis_USA_2_summary.csv')
all_ind = all_ind.sort_values(by=['Pop-celine', 'region-geo', 'taxon', 'ADN-ID'])

pnw_topofile = 'phylonetworks/2020_04_vusa_best_net_ret0.txt'
pnw_topos = import_pnwtopo(pnw_topofile)
# pnw_topos.to_csv('phylonetworks/2020_04_vusa_best_net_ret0.csv')
new_topolist = []
for top in pnw_topos['topo']:
    new_topolist.append(topologic(top, COL_DIC))
pnw_topos['std_topo'] = new_topolist

neww_topolist = []
for top in new_topolist:
    neww_topolist.append(topologic(brootal(top, COL_DIC), COL_DIC))
pnw_topos['rooted_topo'] = neww_topolist

xaxis = all_ind['ADN-ID'].values
all_ind['Pop-celine'].values
BIGDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonet/2020_04_vitis_USA_2_summary_20_reps'


tlist = standardsort(BIGDIR, xaxis, sortlist=None, pnw_topo=pnw_topos['rooted_topo'])
standardsort(BIGDIR, xaxis, sortlist=tlist, pnw_topo=pnw_topos['rooted_topo'])
given_topo = pnw_topos['topo'].values[i]
[rearrange_topo(ptop) for ptop in list(set(neww_topolist))]




tropo = '((vitisusa,vitisasia):0.6301385282681683,(vsylvwest,#H7:0.46259805849744345):0.3100815333675845,(vsylveast,(vveast,(vvwest)#H7:0.5374019415025566):0.163754294985197):0.09807590349937159);'

topologic(brootal(topologic(tropo, COL_DIC), COL_DIC), COL_DIC)


