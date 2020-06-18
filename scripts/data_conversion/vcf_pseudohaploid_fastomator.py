#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

this script enables to transform SNP data in vcf files to fasta format.
it is thought to be used for the Mytilus dataset.


"""
__author__ = 'thomas m huber'
__mail__ = 'thomas.huber@evobio.eu'

#import os
import pandas as pd
import argparse

### use a single vcf field and choose a base due to the given probabilities of the set.
def random_base_set(field,
                    bases):
    import random
    probability = [float(i) for i in field.split(':')[3].split(',')]
    choice = random.random()
    if './.' in field.split(':')[0]:    # sort out uncovered alleles as N
        print('it happened...')
        haplobase = 'N'
    else:
        if choice <= probability[0]:    # choose reference base by random choice
            haplobase = bases[0]
        elif choice <= sum(probability[0:2]):   # choose for heterozygotes
            new_choice = random.random()
            if new_choice <= .5:
                haplobase = bases[0]
            else:
                haplobase = bases[1]
        else:   # choose alternative base if random choice
            haplobase = bases[1]
    return haplobase

### especially for the .
def make_header(atoms):
    pops = ['M.eduEuNorth' if 'M-edu-Europe-N' in i else i for i in atoms]
    pops = ['M.eduEuSouth' if 'M-edu-Europe-S' in i else i for i in pops]
    pops = ['M.eduAm' if 'M-edu-America_' in i else i for i in pops]
    pops = ['M.galloAtl' if 'M-galloprovincialis-A' in i else i for i in pops]
    pops = ['M.galloMed' if 'M-galloprovincialis-M' in i else i for i in pops]
    pops = ['M.tross' if 'M-trossulus' in i else i for i in pops]
    header = [f'>{pop}|{atom}'for pop, atom in zip(pops, atoms)]
    return header


# wordir = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/gitlab/phylogenetwork/data/real_data/mytilus'
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", action="store", dest="vcf_file", required=True,
                    help="Name of the input VCF file")

args = parser.parse_args()
vcf_file = args.vcf_file


print(f'~~~ Conversion of {vcf_file} to fasta format started ~~~')
print(f'~ translation is running in ...~')
with open(vcf_file) as coffee:
    soup = coffee.readlines()   # soup is the list of lines of the file
    sequences = pd.DataFrame([])
    for weight, molecules in enumerate(soup):   # weight reflects the count of the lines
        if weight % 10000 == 0:
            print(f'... line {weight} ...')
        if weight == 15:    # this is the header line of the vcf
            atoms = molecules.strip().split('\t')[9:]
            header = make_header(atoms)
        if weight >= 16:    # these are all the lines with snp data
            bases = molecules.strip().split('\t')[3:5]
            fields = molecules.strip().split('\t')[9:]  # these fields reflect compacted information of one snp of one ind
            basealign = []
            for count, field in enumerate(fields):
                haplobase = random_base_set(field, bases)
                basealign.append(haplobase)
            baseframe = pd.Series(data=basealign)
            baseframe.keys = header
            sequences[weight-16] = baseframe

''.join(sequences.loc[0])
with open(f'{vcf_file.split(".vcf")[0]}.fasta', 'a') as finalfasta:
    for count, head in enumerate(header):
        finalfasta.write(f'{head}\n')
        finalfasta.write(f'{"".join(sequences.loc[count])}\n')

print(f'~~~ Finished: {vcf_file} was successcully converted to fasta format: {vcf_file.split(".vcf")[0]}.fasta. ~~~')


