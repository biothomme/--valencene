#!/usr/bin/env R
# 
# monophorce.r
#
# author = thomas huber
# mail = thomas.huber@evobio.eu
#
# script to infer the branch lengths and topologies for each taxon within a given snp set - we force monophyly!!!
# input: (1) folder, which has SNP fasta files of replicates
#
# output: nexus files of bionj trees of hamming distances for each taxon within dataset
#


library(ape)
library(phangorn)
library(stringr)

##### constants #####
ARGS = commandArgs(trailingOnly = TRUE)

setwd(ARGS[1])
# setwd('/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonet/2020_04_vitis_USA_2_summary_20_reps')
all_dirs = list.dirs()
rep_dirs = all_dirs[grep('rep_', all_dirs)]
for (rep_dir in rep_dirs){
  for (s_file in c('/rep_std_asgnmt.fasta', '/rep_new_asgnmt.fasta')){
    sp_file = paste(rep_dir,s_file, sep = '')
    type_an = str_split(s_file, '_')[[1]][2]
    phyo = read.phyDat(sp_file, format = "fasta", type = "DNA")
    splits = sapply(names(phyo), function(xx) strsplit(xx, '[|]'))
    
    new_labels = c()
    for (spl in splits){
      pasta = strsplit(spl[1], '[.]')
      tax = tolower(paste(pasta[[1]], collapse = ''))
      
      rasta = strsplit(spl[2], '[_]')
      ind = as.character(str_sub(tolower(paste(rasta[[1]], collapse = ''))[1], -10, -1))
      indl = setNames(tax, ind)
      new_labels = c(new_labels, indl)
    }
    names(phyo) = names(new_labels)
    for (taxon in unique(new_labels)) {
      taxsel = which(new_labels == taxon)
      taxset = subset(phyo, subset = taxsel)
      if (taxon == 'vitisusa2'){
        taxon = 'vitisusa'
      }
      file_name = paste(paste(strsplit(sp_file, '/')[[1]][1:2], collapse = '/'), '/topo_', type_an, '_', taxon, '.nex', sep = '')
      
      if (length(taxset) > 1) {
        mat = dist.hamming(taxset, ratio = TRUE, exclude = "none")
        trset = bionj(mat)
        trset = multi2di(trset)
        wrong_way = sapply(trset$tip.label, function(x) paste(rev(strsplit(as.character(x), '')[[1]])[1:9][!is.na(rev(strsplit(as.character(x), '')[[1]])[1:9])], collapse = ''))
        rightway = sapply(wrong_way, function(x) paste(rev(strsplit(x, '')[[1]]), collapse = ''))
        trset$tip.label = unname(rightway)
        write.nexus(trset, file = file_name, translate = F)
      } else {
        wrong_way = sapply(names(taxset), 
                           function(x) paste(rev(strsplit(as.character(x), '')[[1]])[1:9][!is.na(rev(strsplit(as.character(x), '')[[1]])[1:9])], 
                                             collapse = ''))
        rightway = sapply(wrong_way, function(x) paste(rev(strsplit(x, '')[[1]]), collapse = ''))
        
        single_str = paste('[&R] (', unname(rightway), ');', sep = '')
        write(single_str, file = file_name)
      }
    }
  }
}
