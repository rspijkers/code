# This script is meant to be executed from inside the QualityControl environment
# It aims to execute a quality control workflow
# 
# Author: Rik Spijkers (r.spijkers@cern.ch)
# 
# Run this script like this: `alienv setenv QualityControl/latest -c bash qcworkflow.sh`

echo "" # newline

# First things first, let's check if there is a valid token so we can access the ccdb:
alien-token-info > /dev/null # dump output, we are only interested in the exit code
if [ $? -ne 0 ]; then
    while true; do
        read -p "alien-token-info returned with non-zero exit code, do you wish to create a new one? (Y/n)" yn
        case $yn in
            [Yy]* ) alien-token-init; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi

# rm filelist.dat
# # for path in `alien_find -r /alice/data/2022/LHC22o/527228/raw/ [0-9]{4}/o2_ctf.*.root | head -10`
# # for path in `alien_find -r /alice/data/2022/LHC22o/527228/raw/1210/o2_ctf_run00527228_orbit03002.*.root`
# for path in `alien_find -r /alice/data/2022/LHC22o/527228/raw/0850/o2_ctf_run00527228_orbit01699.*.root`
# do
#     echo alien://${path} >> filelist.dat
# done

# make sure we pick the right json file
JSON="json:///home/rik/alice/QualityControl/Modules/ITS/itsTrack.json"
HOST=`hostname`
if [[ "$HOST" == stbc* ]]; then # it seems we are on one of the stbc nodes at nikhef
    echo "we are on stbc"
    JSON="json:///data/alice/rspijker/alice/QualityControl/Modules/ITS/itsTrack.json"
elif [[ "$HOST" == fiora ]]; then
    echo "you are on the Nikhef login server, don't run O2 here!!"
    exit
fi

# let's define some settings here
RAM=8000000000 # 8 gigabyte
# INPUT="903.dat" # can be a filepath, or a txt file with one path per line
INPUT="alien:///alice/data/2023/LHC23zc/537903/raw/1000/o2_ctf_run00537903_orbit0349958816_tf0002069728_epn059.root" # PROBLEMATIC CTF!!! 
# INPUT="alien:///alice/data/2023/LHC23zc/537903/raw/1000/o2_ctf_run00537903_orbit0350641024_tf0002091047_epn180.root" # BEAM SPLASH?
LOG="qcworkflow.log"

echo We are about to execute the workflow, std out is being redirected to $LOG ...
### TRACKS FROM CTF'S ###
o2-ctf-reader-workflow -b --onlyDet ITS --ctf-input $INPUT --remote-regex "^alien:///alice/data/.+" --copy-cmd no-copy --delay 1 --shm-segment-size $RAM | \
o2-its-reco-workflow -b --clusters-from-upstream --disable-mc --trackerCA --tracking-mode async --pipeline its-track:3 --configKeyValues "fastMultConfig.cutMultClusLow=-1;fastMultConfig.cutMultClusHigh=-1;fastMultConfig.cutMultVtxHigh=-1;ITSVertexerParam.phiCut=0.5;ITSVertexerParam.clusterContributorsCut=3;ITSVertexerParam.tanLambdaCut=0.2" | \
o2-qc -b --config $JSON --run > $LOG

### CLUSTERS FROM CTF'S ###
# o2-ctf-reader-workflow -b --onlyDet ITS --ctf-input $INPUT --remote-regex "^alien:///alice/data/.+" --copy-cmd no-copy --delay 1 --shm-segment-size $RAM | \
# o2-qc -b --config $JSON --run > $LOG

# ### Basic workflows to process clusters from TF: ###

# o2-raw-tf-reader-workflow -b --delay 0.2 --input-data 531781_red.dat --max-tf -5 --onlyDet ITS --shm-segment-size 30000000000 | \
# o2-itsmft-stf-decoder-workflow -b --nthreads 5 | \
# o2-qc -b --config json:///home/its/iravasen/QualityControl/Modules/ITS/itsCluster.json --run | tee tracklog.log

# ### Basic workflows to process clusters from CTF: ###

# o2-ctf-reader-workflow -b --onlyDet ITS --ctf-input 529310.dat --remote-regex "^alien:///alice/data/.+" --copy-cmd no-copy --delay 1 --shm-segment-size 30000000000 | \
# o2-qc -b --config json:///home/its/iravasen/QualityControl/Modules/ITS/itsCluster.json --run | tee tracklog.log

# ### Basic workflow to process tracks from TF: ###

# o2-raw-tf-reader-workflow -b --delay 0.2 --input-data 531781_red.dat --max-tf -5 --onlyDet ITS --shm-segment-size 30000000000 | \
# o2-itsmft-stf-decoder-workflow -b --nthreads 5 | \
# o2-its-reco-workflow -b  --clusters-from-upstream --disable-mc --trackerCA --tracking-mode async --configKeyValues "fastMultConfig.cutMultClusLow=-1;fastMultConfig.cutMultClusHigh=-1;fastMultConfig.cutMultVtxHigh=-1;ITSVertexerParam.phiCut=0.5;ITSVertexerParam.clusterContributorsCut=3;ITSVertexerParam.tanLambdaCut=0.2;" | \
# o2-qc -b --config json:///home/its/iravasen/QualityControl/Modules/ITS/itsTrack.json --run | tee tracklog.log
