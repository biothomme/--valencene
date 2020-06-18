#!/usr/bin/env R
# 
# r_bonferroni.R
#
# author = thomas huber
# mail = thomas.huber@evobio.eu
#
# script to cluster populations by their given IBS scores in 4 different ways
# input: (1) folder, which has dist matrix for sat and sylv, as well and loc_mat for sylv, 
#        (2) as well as summary_ind.csv
# output: csv and plots on bionj tree with clusters in color
#

# setwd('/Users/Thomsn/Downloads')

library(ape) 
library(phytools) 
library(gplots) 

##### constants #####
ARGS = commandArgs(trailingOnly = TRUE)
COLS = c('red',
         'yellow', 
         'green', 
         'blue', 
         'orange', 
         'azure4', 
         'bisque3', 
         'palegreen4',
         'turquoise',
         'orchid2',
         'tan4')
SPECIES = c('sativa', 'sylvestris')
METHODS = c('ward.D', 'ward.D2', 'complete', 'mcquitty')
CLUSTERNUM = c(13, 10)

##### functions #####
# obtain geographical color frame
get_color = function(assignm, colors, tree, ind_info, geo = FALSE){
  cuc_mode = rep('', length(assignm))
  for (struc in levels(assignm)){
    lenni = length(cuc_mode[assignm == struc])
    if (geo){
      cuc_mode[assignm == struc] = struc
    } else {
      cuc_mode[assignm == struc] = rep(as.character(colors[struc]), lenni)
    }
  }
  x = as.phylo(tree)
  chaos = x$tip.label
  for (tax in x$tip.label){
    if (length(cuc_mode[ind_info[,2] == tax]) != 1) {
      chaos[chaos == tax] = 'black'
    } else {
      colo = cuc_mode[ind_info[,2] == tax]
      chaos[chaos == tax] = c(colo)
    }
  }
  setNames(chaos, x$tip.label)
}

# plot trees with tree and color
trawdree = function(tree, colors, method_used){
  plot(tree, 
       type = 'r', 
       tip.color = colors,
       show.node.label = TRUE,
       cex = .5,
       lwd=.2,
       align.tip.label = TRUE,
       main = method_used)
}

# get factor for cut to obtain requested amount of clusters for specific hclust
clusterfac_magic = function(hier_clust_result, num_of_clusts){
  fact = 1
  clusts = 0
  while (clusts != num_of_clusts){
    if (clusts > num_of_clusts){
      fact = fact / 1.2
    } else {
      fact = fact * 1.1
    }
    mycl <- cutree(hier_clust_result, h=max(hier_clust_result$height)/fact)
    clusts = length(unique(mycl))  
  }
  fact
}
# obtain vector for coloring the clusters on tree
cluster_color = function(clu_list){
  wdmycol <- rainbow(length(unique(clu_list)), start=0.1, end=0.9) 
  wdmycol[as.vector(clu_list)] 
}

##### main #####
all_clusters = list()
for (sp in SPECIES){
  if (sp == 'sativa'){
    ibfile = paste(c(as.character(ARGS[1]), '/ibs_mat_sativa.csv'), collapse = '')
    datafile = as.character(ARGS[2])
    ibs = as.matrix(read.csv(ibfile, header = TRUE)) #sativa
    ind_data = read.csv(datafile, header = TRUE) #sativa
    struc_mode = setNames(ind_data[,5],ind_data[,2]) #sativa
    clu_num = CLUSTERNUM[1] #sativa
  } else {
    ibfile = paste(c(as.character(ARGS[1]), '/ibs_mat_sylvestris.csv'), collapse = '')
    datafile = paste(c(as.character(ARGS[1]), '/loc_mat_sylvestris.csv'), collapse = '')
    ibs = as.matrix(read.csv(ibfile, header = TRUE)) #sylv
    ind_data = read.csv(datafile, header = TRUE) #sylv
    struc_mode = setNames(ind_data[,3],ind_data[,2]) #sylv
    clu_num = CLUSTERNUM[2] #sylv
  }
  rownames(ibs) = colnames(ibs)[c(2:ncol(ibs))] #sativa
  ibs_new = ibs[,c(2:ncol(ibs))] #sativa
  
  dm = as.dist(ibs_new)
  cols = setNames(COLS[c(1:length(levels(struc_mode)))], levels(struc_mode))

  bn = bionj(dm)
  mush_tree = bn # here the backbone tree can be adjusted - bn is the bionj topology
  pdf(paste(c(as.character(ARGS[1]), sp,'_',clu_num,'_clust.pdf'), collapse = ''))
  geocol = get_color(struc_mode, cols, mush_tree, ind_data)
  geogeo = get_color(struc_mode, cols, mush_tree, ind_data, geo = TRUE)
  clusta_jumbo = data.frame('adn_id' = names(geogeo),
                            'region' = unname(geogeo))
  for (met in METHODS){
    wdc = hclust(dm, method = met)
    wdmycl = cutree(wdc, h=max(wdc$height)/clusterfac_magic(wdc, clu_num)) # 8 is great for sativa!!! 13 for sylv
    wdmycol = cluster_color(wdmycl)
    trawdree(mush_tree, wdmycol, met)
    clusta_jumbo[met] = unname(wdmycl)
  }
  all_clusters[[sp]] = clusta_jumbo  
  
  trawdree(mush_tree, geocol, 'geographic regions')
  add.simmap.legend(colors=cols, prompt=F, fsize=.4, vertical = T, shape = 'circle', y= .75, x = -.125)
  dev.off()
}
# generalize the cluster over all species, to output it
for (each_sp_clust in all_clusters){
  if (exists("clust_info")){
    cur_mac_clust = as.numeric(apply(clust_info, 2, max)[3])
    inc_clust = cbind(each_sp_clust[, c(1:2)],
                      each_sp_clust[, c(3:(2 + length(METHODS)))] + cur_mac_clust)
    clust_info = rbind(clust_info, inc_clust)
  } else{
    clust_info = each_sp_clust
  }
}
# extract the csv including all cluster
output_name = paste(c(as.character(ARGS[1]), 'clusters_vitis_',CLUSTERNUM[1],'sativa_',CLUSTERNUM[2],'sylvestris.csv'), collapse = '')
write.csv(clust_info, 
          file = output_name,
          quote = FALSE,
          row.names = FALSE)

