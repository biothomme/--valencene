#!/bin/bash
# run the script run_america and change name of output
SCRIPT="./scripts/hyde/cluster_run/run_americyde.sh"
JOB=`qsub $SCRIPT | grep -E -o '[0-9]+'`
NUMERIC_ID="./clusterlog/americhyde"
qalter -o ${NUMERIC_ID}.log $JOB
