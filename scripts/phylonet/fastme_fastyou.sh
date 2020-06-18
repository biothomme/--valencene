#!/bin/bash
#
# this short script is made to run FastME on given 
# topologies and alignments for many reps
# input: single rep directory with alignment and 
#        topology

FASTME='../../../softwares/FastME/src/fastme'
TYPE_AN='new'

INFILE=$1/'rep_'$TYPE_AN'_asgnmt.phy' 
TOPO=$1/'forced_mp_wobl_'$TYPE_AN'.newick' 

TREE=$1/'tree_topo_'$TYPE_AN'.new' 
MAT=$1/'tree_mat_'$TYPE_AN'.txt' 
INFO=$1/'tree_info_'$TYPE_AN'.txt' 

$FASTME -i $INFILE -u $TOPO -o $TREE -O $MAT -I $INFO -d 4 -w B
echo '~~~~~ FastME of '$1' is done ~~~~~'

