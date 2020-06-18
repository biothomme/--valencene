#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:05:40 2020

@author: Thomsn
"""
import pandas as pd
challenge = '(vveast,(vvwest,((vsylvwest,tree),((vsylveast,((vitisusa,oooo),vitisasia):0.70920731083779):0.0)#H7:0.03624175498373081::0.8128211905665311):0.20499348795023922):1.9599470411056783,#H7:0.0::0.18717880943346887);'

def rm_nums(topo):
    import re
    nums = r'::?[0-9]\.[0-9]*'
    topo = re.sub(nums, '', topo)
    return topo

def count_str(lst, strng):
    return [el.count(strng) for el in lst]

def split_hyb_start(topo):
    topo = topo.replace(',#H7', ',start')
    return topo

def split_hyb_end(topo):
    work_str = f'{topo.split(")#H7")[0]})'
    postfix = topo.split(")#H7")[1]
    work_str = work_str.split(',')
    op = count_str(work_str, ')')
    cl = count_str(work_str, '(')
    judge = 0
    for i in range(len(op)-1, -1, -1):
        judge = judge + op[i] - cl[i]
        if judge <=0:
            end = i
            break
    prefix = ','.join(work_str[:end])
    main = f'({",".join(work_str[end:])},end)'
    if prefix != '':
        main = f',{main}'
    topo = prefix + main + postfix
    return topo


def remove_too_many_closers(topo):
    cl_deficit = topo.count('(') - topo.count(')')
    if cl_deficit < 0:
        topo = ')'.join(topo.split(')')[:cl_deficit])
    return topo

def remove_doublebr(topo):
    topo = list(topo)
    powerplug = True
    while powerplug:
        if topo[0] == '(' and topo[1] == '(':
            if topo[-1] == ')' and topo[-2] == ')':
                topo = topo[1:-1]
            else:
                powerplug = False
        else:
            powerplug = False
    return ''.join(topo)

def remove_mlessbr(topo):
    if ',' not in topo:
        topo = list(topo)
        powerplug = True
        while powerplug:
            if topo[0] == '(' and topo[-1] == ')':
                    topo = topo[1:-1]
            else:
                powerplug = False
        topo = ''.join(topo)
    return topo

def split_topo(topo):
    work_str = topo.split(',')
    op = count_str(work_str, ')')
    cl = count_str(work_str, '(')
    judge = 0
    end = []
    el = []
    for i in range(len(work_str)):
        judge = judge + cl[i] - op[i]
        el.append(work_str[i])
        if judge <= 0:
            el = ','.join(el).replace(';','')
            el = remove_too_many_closers(el)
            el = remove_mlessbr(el)
            end.append(el)
            el = []
        if judge < 0:
            break
    # topo = ','.join(work_str[:end+1])
    return end

def eat_topo(topo):
    work_str = topo.split('(')
    lauch = []
    for i in range(len(work_str)):
        spl_str = '('.join(work_str[i:])
        lauch.append(split_topo(spl_str))
    return lauch

def all_tax(topo):
    import re
    nam = r'[a-z]*'
    tax_list = [tax for tax in re.findall(nam, topo) if tax != '']
    return tax_list

def topo_framer(topo):
    tax_list = all_tax(topo)
    split_list = eat_topo(tree)
    topo_frame = pd.DataFrame(columns = ['topos', 'ind', 'branches'])
    for node in split_list:
        if len(node) > 1:
            branch_tax = []
            nind = 0
            for branch in node:
                taxa = [tax for tax in tax_list if tax in branch]
                branch_tax.append(taxa)
                nind += len(taxa)
            topo_frame.loc[len(topo_frame)] = [branch_tax, nind, len(branch_tax)]
    return topo_frame

def subset_list(lst):
    return [i for sublist in lst for i in sublist]

def check_all_li_in_li(all_li, in_li):
    return all([True if el in in_li else False for el in all_li])

def check_any_li_in_li(all_li, in_li):
    return any([True if el in in_li else False for el in all_li])

def check_brc_comma_cnt(topo):
    co = topo.count(',')
    op = topo.count('(')
    cl = topo.count(')')
    if (co - cl <= 0 and co - cl <= 0):
        return True
    else:
        return False

def og_root(topoframe, outgroup):
    topo_frame = topoframe.sort_values(by = 'ind')
    og_string = False
    for _, row in topo_frame.iterrows():
        leaves = subset_list(row['topos'])
        if check_all_li_in_li(outgroup, leaves):
            leaves = sorted(row['topos'])
            lvs = []
            for leaf in leaves:
                if check_any_li_in_li(leaf, outgroup):
                    if len(leaf) > 1:
                        lf = f'({",".join(leaf)})'
                    else:
                        lf = leaf[0]
                    lvs.append(lf)
            og_string = f'({",".join(lvs)})'
            break
    return og_string

def update_frame(topo_frame, outgroup):
    for i, row in enumerate(topo_frame['topos']):
        nrow = []
        for li in row:
            lo = [e for e in li if e not in outgroup[1:]]
            lo = [e for e in lo if e not in ['ogsimpsn']]
            lo = ['ogsimpsn' if e == outgroup[0] else e for e in lo]
            if len(lo) > 0:
                nrow.append(lo)
        if (subset_list(nrow) == ['ogsimpsn'] or subset_list(nrow) == ['']):
            nrow = ''
        tree_frame['topos'][i] = nrow
    return tree_frame[tree_frame['topos'] != ''].copy()

def get_ogs(topo_frame):
    og_bool = []
    for row in topo_frame['topos']:
        og = True
        for el in row:
            if 'ogsimpsn' in el:
                og = False
                break
        og_bool.append(og)
    topo_frame['ind'] = [len(subset_list(x)) for x in topo_frame['topos']]
    topo_frame = topo_frame.copy().loc[og_bool,:].sort_values(by = 'ind', ascending = True)
    return topo_frame

def get_ins(topo_frame, ing):
    og_bool = []
    for row in topo_frame['topos']:
        row_list = subset_list(row)
        og_bool.append(check_all_li_in_li(ing, row_list))
    topo_frame['ind'] = [len(subset_list(x)) for x in topo_frame['topos']]
    topo_frame = topo_frame.copy().loc[og_bool,:]
    return topo_frame.sort_values(by = 'ind', ascending = True)

def sort_row(row_list):
    return sorted([sorted(x) for x in row_list])
            
#with open('splitnet1.txt', 'w') as sn1:
with open('outgroups.csv', 'w') as lf1:
    for co, challenge in enumerate(net1_frame['topo']):
        tree = challenge
        tree = rm_nums(tree)
        tree = split_hyb_start(tree)
        tree = split_hyb_end(tree)
#        sn1.write(f'{tree}\n')
        tree_frame = topo_framer(tree)
        outg = ['vitisusa', 'vitisasia']
        og_string = og_root(tree_frame, outg)
        outg = all_tax(og_string)
        if len(outg) < 3:
            lf1.write(f'{",".join(outg)},\n')
        else:
            lf1.write(f'{",".join(outg)}\n')
        
        tree_frame = update_frame(tree_frame, outg)
        ins_frame = get_ogs(tree_frame)
        for row in ins_frame:
            


