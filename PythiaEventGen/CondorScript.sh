#!/bin/bash

# This script is passed to the HTCondor system on Nikhef. It takes an integer as JobID/runnr

source /cvmfs/alice.cern.ch/etc/login.sh
cd /project/alice/users/rspijker/code/PythiaEventGen # Nikhef

ENVIRONMENT="VO_ALICE@pythia::v8304-23,VO_ALICE@ROOT::v6-26-10-alice5-2" # Nikhef

OUTPATH=$1
PYTHIACONFIG=$2
ID=$3

# make outputdir before running script
mkdir -p $OUTPATH

if $? -ne 0; then
  echo "Error: Something went wrong while creating output directory $OUTPATH. Please check the path and permissions."
  exit 1
fi

# run the script, sending the output to tmpdir
eval "alienv setenv $ENVIRONMENT -c ./cascade_correlations $TMPDIR/$ID.root $PYTHIACONFIG $ID"

# copy the output to the final destination
mv $TMPDIR/$ID.root $OUTPATH/$ID.root

# warn the user in case something went wrong
if $? -ne 0; then
  echo "Error: Something went wrong with the script. Please check the output for more details."
fi