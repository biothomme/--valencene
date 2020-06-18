#!/bin/bash
#
# this script starts n jobs for the n different sampling replicates
# of SNPs2CF phylonetwork runs...
#
# input needed: bigdir with repdirs

BIGDIR=$1
EXDIR=`ls $BIGDIR | grep -E rep_`
CURDIR=`pwd`

for DIR in $EXDIR
do
    # here I could input an if statement to not repeat the run
	SCRIPT="${CURDIR}/scripts/phylonetworks/run_july_run.sh"
	JOB=`qsub -b y -N pwn_samrep -cwd -l s_rt=100:00:00 -l h_rt=100:01:00 -j y -pe mpi 10 $SCRIPT "$BIGDIR/$DIR" | grep -E -o '[0-9]+'`
    NUMERIC_ID="${CURDIR}/clusterlog/pnw_new_40"
	qalter -o ${NUMERIC_ID}${DIR}.log -j y $JOB
done
