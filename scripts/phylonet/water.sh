#!/bin/bash
#
# author: thomas markus huber
# contact: thomas.huber@evobio.eu / thomasmarkus.huber.3696@student.uu.se
#
# this script is made to make replicates of phylonet nexus files
# with different seeds using water_the_seeds.py.
# input is a big directory with different replicate directories,
# which have:
#       rep_std_data_{num}_ret.nexus
#       rep_new_data_{num}_ret.nexus
# output are multiple seed replicate files with increasing numbers
# arg1: directory
# arg2: int for num of replicates

BIGDIR=$1
SUBDIRS=`ls $BIGDIR`

WATERSEEDS='./scripts/phylonet/water_the_seeds.py'

FILE1='rep_std_data_0_ret.nexus'
FILE2='rep_new_data_0_ret.nexus'
TOPO1='scen_class_bl_std.new'
TOPO2='scen_class_bl_new.new'


for DIR in $SUBDIRS
do 
	for i in $(seq 1 1 $2)
	do
		python3 $WATERSEEDS -i $BIGDIR$DIR/$FILE1 -t $BIGDIR$DIR/$TOPO1
		python3 $WATERSEEDS -i $BIGDIR$DIR/$FILE2 -t $BIGDIR$DIR/$TOPO2
	done
done
