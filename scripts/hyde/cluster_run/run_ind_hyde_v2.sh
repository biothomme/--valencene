#!/bin/bash
# run the script run_america and change name of output
SCRIPT="./scripts/hyde/cluster_run/ind_hyde_vitis_2.sh"
JOB=`qsub $SCRIPT | grep -E -o '[0-9]+'`
NUMERIC_ID="./clusterlog/runamerichyde"
qalter -o ${NUMERIC_ID}.log $JOB
