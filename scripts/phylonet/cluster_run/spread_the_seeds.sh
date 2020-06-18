#!/bin/bash
#
# this script starts n jobs for the n different seed replicates
# of phylonet runs...
#
# input needed: bigdir with repdirs

BIGDIR=$1
EXDIR=`ls $BIGDIR | head -1`
CURDIR=`pwd`

THATDIR=`ls $BIGDIR$EXDIR'/seed_rep_new_mnr1_ret0' | grep -E ret_ | grep -E nexus`
for FILE in $THATDIR
do
	SCRIPT="${CURDIR}/scripts/phylonet/cluster_run/get_water_flowing.sh"
	JOB=`qsub -b y -N vitis_runpn_rep -cwd -l s_rt=99:00:00 -l h_rt=99:01:00 -j y -pe mpi 10 $SCRIPT $BIGDIR "seed_rep_new_mnr1_ret0/$FILE" | grep -E -o '[0-9]+'`
    NUMERIC_ID="${CURDIR}/clusterlog/pn_mlbm_2"
	qalter -o ${NUMERIC_ID}${FILE}.log -j y $JOB
done

