#!/bin/bash
# this script can be used to start treemix crazyruns for all
# replicates, that have the files within the dir

TRMX="scripts/treemix/crazy_runs.sh"
for DIR in `ls $1`
do
    WDIR="data/real_data/vitis/treemix/permutations/${DIR}/vitis_treemix.txt"
    if test -f $WDIR
    then 
        JOB=`qsub -b y -N treemax -cwd -l s_rt=100:00:00 -l h_rt=100:01:00 -j y $TRMX $WDIR | grep -E -o '[0-9]+'`
        NUMERIC_ID="${CURDIR}/clusterlog/trmax_std_"
        qalter -o ${NUMERIC_ID}${DIR}.log -j y $JOB
        echo $DIR 
    fi 
done

