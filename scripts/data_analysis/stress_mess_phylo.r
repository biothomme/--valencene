setwd('/Users/Thomsn/Desktop/Studium/MEME Programme [current]/Universite de Montpellier [Sem2]/internship_scornavacca_lab/git_lab/phylogenetto/data/real_data/vitis/phylonetworks')

filepath = 'splitnet1.txt'
ogs = read.csv('outgroups.csv', header = F, sep = ',')
library(ape) 
library(phytools) 
library(gplots) 

con = file(filepath, "r")
i = 0
while ( TRUE ) {
  i = i+1
  line = readLines(con, n = 1)
  tree = read.tree(text = line)
  outgroup = c('vitisasia', 'vitisusa')
  if (i == 8){
    outgroup = c('vitisasia', 'vitisusa', 'start')
  }
  ntree = root.phylo(tree, outgroup = outgroup)
  write.tree(ntree, file = 'newtrees.txt', append = TRUE)
  print(outgroup)
  if ( length(line) == 0 ) {
    break
  }
}

  