#!/usr/bin/env R
# 
# monophorce.r
#
# author = thomas huber
# mail = thomas.huber@evobio.eu
#
# script to infer the branch lengths of each taxon within a given snp set monophyly was forced!
# input: (1) folder, which has SNP fasta files of replicates with may_the_monophyly_be_with_you.ph output
#
# output: newick string of scenario classique with updated branchlengths
#
library(ape)
library(stringr)

# install.packages('ape')
library(ape)
## constants ##

TYPESAN = c('new', 'std')
ARGS = commandArgs(trailingOnly = TRUE)

## process ##
    # the_dir = '/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonet/2020_04_vitis_USA_2_summary_20_reps'
    # trial = read.tree(text = '((((((((((b00er3n:0.16969711,vbg90:0.20550045):0.00127897,(b00er4t:0.03106404,vbg169:0.03323286):0.03266938):-0.01670746,((vbg184:0.04505574,b00f6pt:0.04561887):0.00541217,vbg167:0.04710723):0.02068460):0.00549479,b00er4f:0.09088936):0.11257903,((((b00er8a:0.26751490,b00er56:0.22380821):0.00875600,(b00erv4:0.14150547,b00er2l:0.14832196):0.08287549):-0.02370593,b00erqb:0.20687388):0.04415384,(((b00f6nn:0.24198063,b00erd7:0.24327728):0.04043237,b00eqv7:0.27636797):-0.01636502,(b00er6n:0.18325674,b00errm:0.20355564):0.08328901):0.00837423):0.00533875):0.01060458,(((b00eqzg:0.22022059,(b00er7x:0.19996021,b00er0k:0.18425245):0.02467447):0.00856957,b00er7l:0.20274913):0.04678090,b00erb0:0.21860680):0.00929019):0.03153231,vbg146:0.19943227):-0.02008886,vbg027b:0.19801170):-0.00518696,vbg189:0.23736168):0.01821590,vbg271:0.26075542,vbg143d:0.19789415);')
the_dir = ARGS[1]
setwd(the_dir)

all_dirs = list.dirs()
rep_dirs = all_dirs[grep('rep_', all_dirs)]

for (rep_dir in rep_dirs) {
  s_file = '/replicate_data.csv'
  sp_file = paste(rep_dir, s_file, sep = '')
  ind_data = read.csv(sp_file, header = T, sep = ',', dec = '.')
  for (type_an in TYPESAN){
    treefile = paste(rep_dir, '/tree_topo_', type_an, '.new', sep = '')
    trial = read.tree(file = treefile)
    ### attention: for the phylip-file I cut the filenames to 9 instead of 10 chars! this needs to be accounted here!
    ##### root it #####
    if (type_an == 'std') {
      ind_data$pop = sapply(ind_data$Pop.celine, function(x) gsub('[. 2]', '', str_to_lower(x)))
    }
    subsind = ind_data[ind_data$pop == 'vitisasia' | ind_data$pop == 'vitisusa',]
    wrong_way = sapply(subsind$cut_names, 
                       function(x) paste(rev(strsplit(as.character(x), '')[[1]])[1:9][!is.na(rev(strsplit(as.character(x), '')[[1]])[1:9])], 
                                         collapse = ''))
    rightway = sapply(wrong_way, function(x) paste(rev(strsplit(x, '')[[1]]), collapse = ''))
    subsind$cut_names = unname(rightway)
    for (st in subtrees(trial)){
      if (paste(sort(subsind$cut_names), collapse = '') == paste(sort(st$tip.label), collapse = '')){
        rootnode = st$node.label[1]
      }
    }
    # rootrial = root(trial, outgroup = subsind$cut_names, resolve.root = T)
    rootrial = root(trial, node = rootnode, resolve.root = T)
    ##### get subtrees and taxon brlengths ####
    tax_lens = data.frame(row.names = c('tax', 'branch_length'))
    for (tax in unique(ind_data$pop)){
      subsind = ind_data[ind_data$pop == tax,]
      wrong_way = sapply(subsind$cut_names, 
                         function(x) paste(rev(strsplit(as.character(x), '')[[1]])[1:9][!is.na(rev(strsplit(as.character(x), '')[[1]])[1:9])], 
                                           collapse = ''))
      rightway = sapply(wrong_way, function(x) paste(rev(strsplit(x, '')[[1]]), collapse = ''))
      subsind$cut_names = unname(rightway)
      for (st in subtrees(rootrial)){
        if (paste(sort(subsind$cut_names), collapse = '') == paste(sort(st$tip.label), collapse = '')){
           kt = st
           brlens = data.frame(row.names = c('ind', 'branch_length'))
           for (ed in c(1:length(kt$tip.label))) {
             edsum = 0
             edbool = (kt$edge[,2] == ed)
             
             while (any(edbool)) {
               edval = kt$edge.length[edbool]
               edsum = edsum + edval
               ned = kt$edge[edbool,1]
               edbool = kt$edge[,2] == ned
             }
             brlens[as.character(ed)] = c(kt$tip.label[ed], edsum)
             inlen = mean(sapply(brlens[2,], as.numeric))
             outlen = rootrial$edge.length[rootrial$edge[,2] == st$name]
             taxlen = outlen + inlen
           }
        }
      }
      if (!(tax %in% colnames(tax_lens))) {
        single_node = which(rootrial$tip.label == subsind$cut_names[1])
        sinbool = rootrial$edge[,2] == single_node
        taxlen = rootrial$edge.length[sinbool]
      }
      if (taxlen <= 0) {
        taxlen = '0.00001'
      }
      tax_lens[tax] = c(tax, taxlen)
    }
    
    ##### get intertax branchlengths #####
    VITISCOMBO = c('vitis', 'east', 'west', 'st')
    for (vcb in VITISCOMBO){
      subsind = ind_data[sapply(ind_data$pop, function(x) grepl(vcb, x)),]
      wrong_way = sapply(subsind$cut_names, 
                         function(x) paste(rev(strsplit(as.character(x), '')[[1]])[1:9][!is.na(rev(strsplit(as.character(x), '')[[1]])[1:9])], 
                                           collapse = ''))
      rightway = sapply(wrong_way, function(x) paste(rev(strsplit(x, '')[[1]]), collapse = ''))
      subsind$cut_names = unname(rightway)
      for (st in subtrees(rootrial)){
        if (paste(sort(subsind$cut_names), collapse = '') == paste(sort(st$tip.label), collapse = '')){
          medlen = rootrial$edge.length[rootrial$edge[,2] == st$name]
          if (medlen <= 0) {
            medlen = '0.00001'
          }
          tax_lens[vcb] = c(vcb, medlen)
        }
      }
    }
    
    ##### build scenario classique #####
    pax_lens = tax_lens[2,]
    newick_str = paste('((vitisasia:', 
          pax_lens$vitisasia,
          ',vitisusa:',
          pax_lens$vitisusa,
          '):',
          pax_lens$vitis,
          ',((vsylveast:',
          pax_lens$vsylveast,
          ',vveast:',
          pax_lens$vveast,
          '):',
          pax_lens$east,
          ',(vsylvwest:',
          pax_lens$vsylvwest,
          ',vvwest:',
          pax_lens$vvwest,
          '):',
          pax_lens$west,
          '):',
          pax_lens$st,
          ');',
          sep = '',
          collapse = '')
    out_file = paste(rep_dir, '/scen_class_bl_', type_an, '.new', sep = '')
    write(newick_str, out_file)
  }
}

