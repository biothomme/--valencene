#!/bin/bash
#
# because my folders are becoming messier,
# I need to clean up. Let me start with 
# putting all 100 seed_reps in specific folders.
# Input: bigdir dith repdirs with seedreps


BIGDIR=$1
DIRRO=`ls $BIGDIR | grep -E rep`

echo $DIRRO
EXT1='new'
EXT2='std'

for DIR in $DIRRO
do
	# NEWDIR1=${BIGDIR}${DIR}'/seed_rep_'${EXT1}'_mnr1_ret0'
	# mkdir $NEWDIR1
	# mv ${BIGDIR}${DIR}/*${EXT1}_data_0_ret_*.nexus $NEWDIR1
        # NEWDIR2=${BIGDIR}${DIR}'/seed_rep_'${EXT2}'_mnr1_ret0'
        # mkdir $NEWDIR2
        # mv ${BIGDIR}${DIR}/*${EXT2}_data_0_ret_*.nexus $NEWDIR2

        TOPODIR1=${BIGDIR}${DIR}'/topologies_'${EXT1}'_ret0'
        mkdir $TOPODIR1
        mv ${BIGDIR}${DIR}/topo_${EXT1}* $TOPODIR1
	mv ${BIGDIR}${DIR}/tree_*_${EXT1}* $TOPODIR1
	mv ${BIGDIR}${DIR}/forced_*${EXT1}.newick $TOPODIR1
	mv ${BIGDIR}${DIR}/scen_class_bl_${EXT1}* $TOPODIR1
	mv ${BIGDIR}${DIR}/*${EXT1}*.phy $TOPODIR1
        TOPODIR2=${BIGDIR}${DIR}'/topologies_'${EXT2}'_ret0'
        mkdir $TOPODIR2
        mv ${BIGDIR}${DIR}/topo_${EXT2}* $TOPODIR2
        mv ${BIGDIR}${DIR}/tree_*_${EXT2}* $TOPODIR2
        mv ${BIGDIR}${DIR}/forced_*${EXT2}.newick $TOPODIR2
        mv ${BIGDIR}${DIR}/scen_class_bl_${EXT2}* $TOPODIR2
        mv ${BIGDIR}${DIR}/*${EXT2}*.phy $TOPODIR2

done
