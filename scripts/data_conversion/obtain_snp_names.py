#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 09:36:38 2020

this script extracts the names of the SNPs as a string - for vitis dataset

@author: Thomsn
"""
import os
import pandas as pd

CURDIR = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis'
os.chdir(CURDIR)

data_file = '2020_04_vitis_USA_2.txt'
snps = pd.DataFrame(columns = ['chromosome', 'locus'])
with open(data_file, 'r') as daf:
    soup = daf.readlines()
    for molecule in soup:
        if 'chr' in molecule:
            separ = molecule.split('\t')
            snp = separ[0]
            chromosome = snp.split('_')[0]
            locus = int(snp.split('_')[-1])
            if 'Un' in chromosome:
                chromosome = int(99)
            else:
                chromosome = int(chromosome.split('chr')[-1])
            snps.loc[len(snps),:] = [chromosome, locus]
snps.to_csv('2020_04_vitis_snps_summary.csv')


