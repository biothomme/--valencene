#!/usr/bin/env R
#
# r_bonferroni.R
#
# author = thomas huber
# mail = thomas.huber@evobio.eu
#
# script to analyze HyDe (python package phyde) outputs
# previously set p_values will be corrected after bonferroni
#

args = commandArgs(trailingOnly = TRUE)
ip_name = as.character(args[1L])
dato = read.table(ip_name, header = T)
plato = dato
plato$Pvalue = p.adjust(dato$Pvalue, method = 'bonferroni')
op_name = paste(gsub('.txt','',ip_name),'-bonf.txt', sep = '')
write.table(plato, op_name, sep="\t", quote = FALSE, row.names = FALSE)

