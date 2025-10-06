# This script is meant to be executed from inside the O2 environment
# It aims to execute an O2 workflow, redirecting the std::out to log_o2.txt
# 
# Author: Rik Spijkers (rik.spijkers@cern.ch)
# 
# Run this script like this: `alienv setenv O2Physics/latest -c bash runtest.sh`

# TODO: add option to pass txtfile with AO2D's to run over? maybe pass flags like '--mc' or '--run2' to automatically select the correct filelist?

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

FOLDER="local"
HOST=`hostname`
if [[ "$HOST" == stbc* ]]; then # it seems we are on one of the stbc nodes at nikhef
    echo "we are on stbc"
    FOLDER="stbc"
elif [[ "$HOST" == fiora ]]; then
    echo "you are on the Nikhef login server, don't run O2 here!!"
    exit
fi

# This is obsolete: we normally just test on 22o_apass4 sample
# if we are not on stbc or fiora, assume we are running locally
FILELIST=$FOLDER/LHC22o_apass4_5.txt

# see if we can run over the 22o_apass4 sample

# define the config so we can conveniently pass it to each workflow
CONFIG="-b --configuration json://ssbarconfigMCclosure.json" # -b means no debug gui

echo "hey it seems we are about to execute the workflow, cool!"
# execute the entire workflow
o2-analysis-mccollision-converter $CONFIG | o2-analysis-tracks-extra-v002-converter $CONFIG | \
o2-analysis-timestamp $CONFIG | o2-analysis-track-propagation $CONFIG | o2-analysis-event-selection $CONFIG | o2-analysis-multiplicity-table $CONFIG | \
o2-analysis-pid-tpc-base $CONFIG | o2-analysis-pid-tpc $CONFIG | \
o2-analysis-lf-strangenessbuilder $CONFIG | \
o2-analysis-lf-cascadecorrelations $CONFIG --aod-file "alien:///alice/sim/2024/LHC24j1/526641/AOD/001/AO2D.root" > log_o2.txt

# --resources-monitoring 10
# in case of skimmed datasets: o2-analysis-lf-strangeness-filter
# sometimes you may need this workflow: o2-analysis-collision-converter $CONFIG | 
# when running on older stuff maybe add this: o2-analysis-zdc-converter $CONFIG | 
# when running on older data use (before September 2023?): o2-analysis-bc-converter $CONFIG | 
# when running on older data use (before November 2023?): o2-analysis-tracks-extra-converter $CONFIG | 
# when running on older data use (before December 2023?): o2-analysis-v0converter $CONFIG | 
# @FILELIST
# @24j1_10files.txt

### 22o pass6: alien:///alice/data/2022/LHC22o/526641/apass6/0650/o2_ctf_run00526641_orbit0215379968_tf0000066791_epn157/001/AO2D.root
### 22o pass7?: alien:///alice/data/2022/LHC22o/526641/apass7/0650/o2_ctf_run00526641_orbit0215950848_tf0000071251_epn246/001/AO2D.root
### 22o pass7 skimmed: alien:///alice/data/2022/LHC22o/528543/apass7_skimmed/0950/o2_ctf_run00528543_orbit0590742144_tf0000112425_b7s02p9411/001/AO2D.root
### 24am skimmed: alien:///alice/data/2024/LHC24am/555976/apass1_skimmed/0610/o2_ctf_run00555976_orbit0465019808_tf0000243738_b9p10p3582/001/AO2D.root
### MC file? /alice/sim/2023/LHC23k2c/1/529662/001/AO2D.root
### MC genpurp LHC24b1b (22o_apass6 anchor?): /alice/sim/2024/LHC24b1b/0/528531/AOD/001/AO2D.root

# check the exit code, if it's non-zero then let the user know.
if [ $? -ne 0 ]; then
    echo "The O2 workflow exited with non-zero exit code, something probably went wrong! Please check the log: log_o2.txt."
fi

# let the user know the script is finished
# notify-send "Your O2 analysis is done!" # pop-up DOESNT WORK ON STBC
# echo -e "\a" # sound notification
