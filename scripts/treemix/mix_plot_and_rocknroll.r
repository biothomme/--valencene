#!/usr/bin/env R
#
# mix_plot_and_rocknroll.r
#
# author = thomas huber
# mail = thomas.huber@evobio.eu
#
# script to plot treemix output
#

source("scripts/treemix/plotting_funcs.R")

args = commandArgs(trailingOnly = TRUE)
ip_tree = as.character(args[1L])

op_tree = gsub(".gz", ".pdf", ip_tree)

pdf(op_tree)
plot_tree(ip_tree)
dev.off()

