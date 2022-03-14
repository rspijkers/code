#!/bin/bash

# This script is passed to the batch queue on Nikhef with the qsub command, it aims to run a pythia script

source /cvmfs/alice.cern.ch/etc/login.sh
cd /user/rspijker/code # Nikhef

ENVIRONMENT="VO_ALICE@AliPhysics::vAN-20180913-1,VO_ALICE@pythia::v8230-1" # Nikhef
# ENVIRONMENT="O2/latest,pythia/latest" # Home

# TODO: add outputfile location as option of cpp script? 
# so we can steer it to /dcache/ at nikhef, or something else.
# we can then also let the output path depend on where we run it, home or nikhef
eval "alienv setenv $ENVIRONMENT -c ./ssbar_correlations"
