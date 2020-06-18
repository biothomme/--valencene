#!/bin/bash
# run the script run_america and change name of output
SCRIPT="./scripts/hyde/cluster_run/run_dioli.sh"
JOB=`qsub $SCRIPT | grep -E -o '[0-9]+'`
NUMERIC_ID="./clusterlog/diolo_hyde_am"
qalter -o ${NUMERIC_ID}.log $JOB
