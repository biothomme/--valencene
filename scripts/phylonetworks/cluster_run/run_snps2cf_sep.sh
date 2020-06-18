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
	SCRIPT="${CURDIR}/scripts/phylonetworks/cluster_run/snp2cf_rep_vitis.sh"
	JOB=`qsub -b y -N snpstocf_samrep -cwd -l s_rt=200:00:00 -l h_rt=200:01:00 -j y -pe mpi 10 $SCRIPT "$BIGDIR/$DIR" | grep -E -o '[0-9]+'`
    NUMERIC_ID="${CURDIR}/clusterlog/snps2cf_run2_"
	qalter -o ${NUMERIC_ID}${DIR}.log -j y $JOB
done

