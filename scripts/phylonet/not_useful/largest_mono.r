#!/usr/bin/env R
# 
# largest_mono.r
#
# author = thomas huber
# mail = thomas.huber@evobio.eu
#
# script to find the largest monophylies for each taxon within a given SNP set
# input: (1) replicate folder, which has topology and replicate info files
#
# output: lists of largest monophyletic sets per taxon
#

# MORE OR LESS USELESS for our analysis - we do not have monophyletic taxa...

library(ape)
library(phangorn)
library(stringr)

##### FUNCTIONS #####
### get biggest monophylum of specified group of ind ###
selected_ind = selnam
tree = tr
taxa = direction
catch_monos = function(selected_ind, tree, taxa) {
  # wall_monop = data.frame()
  seltax = names(selected_ind)
  
  monop = data.frame()
  for (i in c(length(selected_ind):1)){
    psblts = combn(seltax, i)
    for (psb in c(1:ncol(psblts))) {
      sel_tax = psblts[,psb]
      check2 = is.monophyletic(tree, sel_tax)
      if (check2) {
        the_line = paste(sel_tax, collapse = ',')
        new_line = data.frame('taxon' = taxa, 'monophylum' = the_line)
        monop = rbind(monop, new_line)
      }
    }
    if (nrow(monop) > 0) {
    #   for (columca in c(1:ncol(monop))){
    #     columna = monop[,columca]
    #     the_line = paste(columna, collapse = ',')
    #     new_line = data.frame('taxon' = taxa, 'monophylum' = the_line)
    #     wall_monop = rbind(wall_monop, new_line)
    #   }
    #   new_line = data.frame('taxon' = taxa, 'monophylum' = the_line)
       break
    }
  }
  monop
}
### for general mps ###
simple_monophylies = function(tree, assigned_names){
  all_monop = data.frame()
  the_tips = assigned_names[tree$tip.label]
  cols = as.numeric(as.factor(the_tips))
  
  for (taxc in unique(cols)) {
    tax = levels(as.factor(assigned_names[tree$tip.label]))[taxc]
    selnam = the_tips[unname(the_tips) == tax]
    all_monop = rbind(all_monop, catch_monos(selnam, tree, tax))
  }
  all_monop
}
### for monophyles with west, east distribution ###
sclas_monophylies = function(tree, assigned_names){
  all_monop = data.frame()
  the_tips = assigned_names[tree$tip.label]
  cols = as.numeric(as.factor(the_tips))
  
  for (direction in c('west', 'east')) {
    meck = sapply(the_tips, function(xx) grepl(direction, xx, ignore.case = T))
    selnam = as.factor(assigned_names[tree$tip.label])[meck]
    all_monop = rbind(all_monop, catch_monos(selnam, tr, direction))
  }
  all_monop
}

##### ANALYSIS #####
ARGS = commandArgs(trailingOnly = TRUE)
setwd(ARGS[1])
# setwd('/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/')

tr = read.nexus('rep_new_topo.nex')
summatr = read.csv('replicate_data.csv')
summatr$cut_names[summatr$cut_names %in% tr$tip.label]
plot(tr)
sumnam = setNames(summatr$pop, summatr$cut_names)

tre_tips = sumnam[tr$tip.label]

colr = as.numeric(as.factor(tre_tips))
plot(tr, tip.color = colr)
write.nexus(tr, file = '')
tr$edge.length

simple_monophylies(tr, sumnam)

sclas_monophylies(tr, sumnam)
