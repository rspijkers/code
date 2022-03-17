#!/bin/bash

# This script steers a pythia simulation on stbc

ENVIRONMENT="VO_ALICE@AliPhysics::vAN-20180913-1,VO_ALICE@pythia::v8230-1" # Nikhef
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
