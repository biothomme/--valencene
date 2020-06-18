#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 08:36:54 2020

this script should convert a nexus file to different xml files, needed for 
snappnet

@author: Thomsn
"""

__author__ = 'thomas m. huber'
__email__ = ['thomas.huber@evobio.eu']

import math
import os
import pandas as pd

def wprint(x, file):
    file.write(x)
    return

def xmlfiler(mcmc_file,
             mcmc_txmplate,
             rep_data,
             nex_file):
    with open(mcmc_file, 'w') as nwf:
        with open(mcmc_txmplate, 'r') as mct:
            mct_soup = mct.readlines()
            for mctolecule in mct_soup[0:4]:
                wprint(mctolecule, nwf)
            wprint(f'id="nexus_{rep.replace("_","")}"\n', nwf)
            for mctolecule in mct_soup[5:7]:
                wprint(mctolecule, nwf)
                # wprint(f'id="nexus_{rep}"\n')
            go2 = False
            for mctolecule in mct_soup[7:]:
                if '</data>' in mctolecule:
                    with open(nex_file, 'r') as nxs:
                        nxs_soup = nxs.readlines()
                        go = True
                        for nxsolecule in nxs_soup[6:]:
                            if ';End;' in nxsolecule:
                                go = False
                            if go:
                                ind = nxsolecule.split(' ')[0]
                                seq = nxsolecule.split(' ')[-1].replace('\n', '')
                                tax = simplifiy_tax(rep_data['pop'][rep_data['cut_names'] \
                                                                    == ind].values[0])
                                ind_str = f'{"".join([" "]*24)}<sequence id="seq_{ind}_{tax}" taxon="{ind}_{tax}" totalcount="2" value="{seq}"/>\n'
                                wprint(ind_str, nwf)
                    go2 = True
                if go2:
                    if not 'taxonset' in mctolecule:
                        wprint(mctolecule.replace('nexus_rep01',f'nexus_{rep.replace("_","")}'), nwf)
                    else:
                        go2 = False
                        for tax in set(rep_data['pop']):
                            ind_list = rep_data['cut_names'][rep_data['pop'] == tax].values
                            tax = simplifiy_tax(tax)
                            tax_str = mctolecule.replace('vve', tax)
                            wprint(tax_str, nwf)
                            for ind in ind_list:
                                ind_str = f'{"".join([" "]*24)}<taxon id="{ind}_{tax}" spec="Taxon"/>\n'
                                wprint(ind_str, nwf)
                            wprint(f'{"".join([" "]*20)}</taxonset>\n', nwf)
                if '</taxa>' in mctolecule:
                    go2 = True
                    wprint(mctolecule, nwf)
    return

def alltrees(taxa):
    import itertools as it
    length = ['0.25',
                '0.25',
                '0.5',
                '0.25',
                '0.25',
                '0.25',
                '0.5',
                '0.125',
                '0.625',
                '0.125',
                '0.25',
                '1.0']
    topo_list = []
    for midg in it.combinations(taxa[2:],2):
        og = f'({taxa[0]}:t01,{taxa[1]}:t02)S1:t03'
        cherry = [tax for tax in taxa[2:] if tax not in midg]
        chg = f'(({cherry[0]}:t04,{cherry[1]}:t05)S2:t06'
        for tax in midg:
            tox = [to for to in midg if to != tax][0]
            mdg = f',{tax}:t07)S3:t08,{tox}:t09)S4:t10'
            total_str = f'{og},({chg}{mdg})S5:t12'
            for i in range(1,13):
                if i < 10:
                    ii = f'0{i}'
                else:
                    ii = i
                tt = f't{ii}'
                total_str = total_str.replace(tt, length[i-1])
            topo_list.append(f'(({total_str})')
    for midg in taxa[3:]:
        og = f'({taxa[0]}:t01,{taxa[1]}:t02)S1:t03'
        cherry = [tax for tax in taxa[3:] if tax != midg]
        chg1 = f'({cherry[0]}:t04,{cherry[1]}:t05)S2:t06'
        chg2 = f'({midg}:t04,{taxa[2]}:t05)S3:t06'
        total_str = f'{og},({chg1},{chg2})S4:t11)S5:t12'
        for i in range(1,13):
            if i < 10:
                ii = f'0{i}'
            else:
                ii = i
            tt = f't{ii}'
            total_str = total_str.replace(tt, length[i-1])
        topo_list.append(f'(({total_str})')
    return topo_list

def statement(ml_tempstate,
              TAXA,
              rep):
    topologies = pd.DataFrame({'name':[f'0{i}' if i < 10 else f'{i}' for i in range(1,16)],
                          'tree':alltrees(TAXA)})
    for _, topo in topologies.iterrows():
        ml_statefile = f'{rep}/SnappNetVigne/ML/MLVigne_{topo["name"]}.xml.state'
        topology = topo['tree']
        with open(ml_statefile, 'w') as nml:
            with open(ml_tempstate, 'r') as mlt:
                mltsoup = mlt.readlines()
                for mltecule in mltsoup:
                    post = ' </statenode>'
                    if 'coalescenceRate' in mltecule:
                        coal_rates = ' '.join([str(math.exp(1))] * topology.count(':'))
                        pre = f"<statenode id='coalescenceRate'>coalescenceRate[{topology.count(':')} 1] (0.0,10.0): "
                        total_str = f'{pre}{coal_rates}{post}\n'
                    elif 'network:species' in mltecule:
                        pre = "<statenode id='network:species'>"
                        total_str = f'{pre}{topology}{post}\n'
                    else:
                        total_str = mltecule
                    nml.write(total_str)
    return


BIGDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/snappnet'
TAXA = ['vas', 'vus', 'vve', 'vvw', 'vse', 'vsw']

os.chdir(BIGDIR)

rep_list = [rep for rep in os.listdir() if 'rep_' in rep]
simplifiy_tax = lambda x: x.replace('ast','').replace('ylv','').replace('itis',\
    '').replace('ia','').replace('usa','us').replace('est','')

for rep in rep_list:
    nex_file = f'{rep}/rep_new_data_0_ret.nexus'
    rep_file = f'{rep}/replicate_data.csv'
    mcmc_txmplate = f'{rep}/SnappNetVigne/MCMC/MCMCVigne.xml'
    ml_txmplate = f'{rep}/SnappNetVigne/ML/MLVigne.xml'
    ml_tempstate = f'{rep}/SnappNetVigne/ML/MLVigne.xml.state.txt'
    
    mcmc_file = f'{rep}/SnappNetVigne/MCMC/MCMCVigne_new.xml'
    ml_tfile = f'{rep}/SnappNetVigne/ML/MLVigne_new.xml'


    rep_data = pd.read_csv(rep_file)    

# THIS PART IS FOR MCMC    
    xmlfiler(mcmc_file, mcmc_txmplate, rep_data, nex_file)
# THIS PART IS FOR ML
    xmlfiler(ml_tfile, ml_txmplate, rep_data, nex_file)
    statement(ml_tempstate, TAXA, rep)
';\n'.join(alltrees(TAXA))
            