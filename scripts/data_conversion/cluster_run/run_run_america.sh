#!/bin/bash
# run the script run_america and change name of output
SCRIPT="./scripts/data_conversion/cluster_run/run_america.sh"
JOB=`qsub $SCRIPT | grep -E -o '[0-9]+'`
NUMERIC_ID="./clusterlog/123"
qalter -o ${NUMERIC_ID}.log $JOB
