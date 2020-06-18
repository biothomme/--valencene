#!/usr/bin/env R
# 
# longest_branch_alive.R
#
# author = thomas huber
# mail = thomas.huber@evobio.eu
#
# script to infer the branch lengths and topologies for given snp set
# input: (1) folder, which has SNP fasta files of replicates
#
# output: nexus files of bionj trees of hamming distances for dataset
#

# setwd('/Users/Thomsn/Downloads')


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
    x = read.phyDat(sp_file, format = "fasta", type = "DNA")
    mat = dist.hamming(x, ratio = TRUE, exclude = "none")
    tree = bionj(mat)
    
    ## the headers need to be adjusted! ##
    splits = sapply(tree$tip.label, function(e) strsplit(e, '[|]'))
    new_labels = c()
    for (spl in splits){
      #  --- taxa not neccessary for label --- #
      # pasta = strsplit(spl[1], '[.]')
      # tax = tolower(paste(pasta[[1]], collapse = ''))
      
      rasta = strsplit(spl[2], '[_]')
      ind = str_sub(tolower(paste(rasta[[1]], collapse = ''))[1], -10, -1)
      new_labels = c(new_labels, ind)
    }
    tree$tip.label = new_labels
    
    ## save nexus file ##
    nex_file = paste(strsplit(sp_file, 'rep_')[1],'_topo.nex', sep = '')
    write.nexus(tree, file = nex_file, translate = FALSE)
  }
}
