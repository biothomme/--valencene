#!/bin/bash
#
# because my folders are becoming messier,
# I need to clean up. Let me start with 
# putting all 100 seed_reps in specific folders.
# Input: bigdir dith repdirs with seedreps


BIGDIR=$1
DIRRO=`ls $BIGDIR | grep -E rep_`

echo $DIRRO
EXT1='new'

for DIR in $DIRRO
do
	NEWDIR1=${BIGDIR}${DIR}'/seed_rep_'${EXT1}'_mnr1_ret1'
	mkdir $NEWDIR1
	mv ${BIGDIR}${DIR}/*${EXT1}_data_1_ret_*.nexus $NEWDIR1
done
