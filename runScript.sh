#!/bin/bash

# This script steers a pythia simulation on stbc

ENVIRONMENT="VO_ALICE@pythia::v8304-9,VO_ALICE@ROOT::v6-24-06-18" # Nikhef
# ENVIRONMENT="O2/latest,pythia/latest" # Home
echo $ENVIRONMENT

# cd /user/rspijker/code # Nikhef
# if [[ ! -d "output" ]]
# then 
#     mkdir output
# fi
# make pythia script
eval "alienv setenv $ENVIRONMENT -c make ssbar_correlations"
for i in {1..10}
do
    echo "submit number $i"
    qsub submitScript.sh
done
