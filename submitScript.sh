#!/bin/bash

# This script is passed to the batch queue on Nikhef with the qsub command, it aims to run a pythia script

source /cvmfs/alice.cern.ch/etc/login.sh
cd /user/rspijker/code # Nikhef

ENVIRONMENT="VO_ALICE@pythia::v8304-9,VO_ALICE@ROOT::v6-24-06-18" # Nikhef
# ENVIRONMENT="O2/latest,pythia/latest" # Home

# get job id and take only the numeric part (leave out the .burrell.nikhef.nl)
JOB_ID=$(echo $PBS_JOBID | cut -d'.' -f 1)

OUTPUTDIR=$1
OUTPUTFILE="$OUTPUTDIR/$JOB_ID.root"
RUNNR=$2

eval "alienv setenv $ENVIRONMENT -c ./ssbar_correlations $OUTPUTFILE $RUNNR"
