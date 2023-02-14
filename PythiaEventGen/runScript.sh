#!/bin/bash

# This script steers a pythia simulation on stbc
cd /project/alice/users/rspijker/code/PythiaEventGen # Nikhef

ENVIRONMENT="VO_ALICE@pythia::v8304-23,VO_ALICE@ROOT::v6-26-10-alice5-2" # Nikhef
# ENVIRONMENT="O2/latest,pythia/latest" # Home
echo $ENVIRONMENT

# define outputdir, check to see if it exists, if not then mkdir it
OUTPUTPATH="/dcache/alice/rspijker/ModelStudyJan"
OUTPUTDIR="Monash_pp_100M_14TeV"
[ ! -d "$OUTPUTPATH/$OUTPUTDIR" ] && mkdir -p "$OUTPUTPATH/$OUTPUTDIR" && echo "outputdir created"
LOGPATH="$OUTPUTPATH/$OUTPUTDIR/logs"
[ ! -d "$LOGPATH" ] && mkdir -p "$LOGPATH"

# pythia config file
PYTHIACONFIG="pythia_settings/ssbar_monash.cmnd"

# make pythia script
alienv setenv $ENVIRONMENT -c make ssbar_correlations -B
if [ $? -ne 0 ]; then
    echo "Compilation was not succesful, aborting..."
    exit
fi

# submit the jobs
for i in {1..10}
do
    echo "submit number $i"
    qsub -F "$OUTPUTPATH/$OUTPUTDIR $PYTHIACONFIG $i" submitScript.sh -o $LOGPATH -e $LOGPATH
done
