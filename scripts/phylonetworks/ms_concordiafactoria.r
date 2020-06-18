#!/usr/bin/env R
# 
# ms_concordiafactoria.R
#
# author = thomas huber
# mail = thomas.huber@evobio.eu
#
# script to input a SNPs2CF analysis. It can be used
# for SNP data and takes advantage of the package
# SNPs2CF by melisa olave (2020).
# this script enables to analyze multiple replicates
# of a dataset. This is needed, because SNPs2CF is
# limited to a maximal amount of individuals per
# quartet (100 000).

# install.packages("foreach")
# install.packages("doMC")

args = commandArgs(trailingOnly = TRUE)
ip_name = as.character(args[1])
mapname = as.character(args[2])
op_name = sub('data.phy', 'out.csv', ip_name)

pack = 'scripts/phylonetworks/SNPs2CF.r'
source(pack)
output = SNPs2CF(seqMatrix=ip_name,
                 ImapName=mapname,
                 between.sp.only = TRUE,
                 max.SNPs = 100000,
                 bootstrap=FALSE, # NOTE FOR LATER: DO WE WANT BOOTSTRAPS?
                 outputName=op_name,
                 save.progress=FALSE)


