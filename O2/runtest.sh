# This script is meant to be executed from inside the O2 environment, perhaps this should be changed
# It aims to execute an O2 workflow, redirecting the std::out to log_o2.txt
# 
# Author: Rik Spijkers (r.spijkers@cern.ch)
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

# TODO: check if we are on local or stbc by bash command `hostname`, select the correct filelist based on this info
FOLDER="local"
HOST=`hostname`
if [ $HOST == stbc* ]; then # it seems we are on one of the stbc nodes at nikhef
    echo "we are on stbc"
    FOLDER="stbc"
elif [ $HOST == fiora ]; then
    echo "you are on the Nikhef login server, don't run O2 here!!"
    exit
fi
# if we are not on stbc or fiora, assume we are running locally
FILELIST=$FOLDER/LHC22t_apass3.txt

# define the config so we can conveniently pass it to each workflow
CONFIG="-b --configuration json://ssbarconfig.json" # -b means no debug gui

echo "hey it seems we are about to execute the workflow, cool!"
# execute the entire workflow
o2-analysis-timestamp $CONFIG | o2-analysis-track-propagation $CONFIG | o2-analysis-event-selection $CONFIG | o2-analysis-multiplicity-table $CONFIG | o2-analysis-pid-tpc-base $CONFIG | o2-analysis-pid-tpc $CONFIG | o2-analysis-lf-lambdakzerobuilder $CONFIG | o2-analysis-lf-cascadebuilder $CONFIG | o2-analysis-lf-cascadeanalysis $CONFIG | o2-analysis-lf-strangecorrelations $CONFIG --aod-file alien:///alice/data/2022/LHC22t/529552/apass3/0800/o2_ctf_run00529552_orbit0007634688_tf0000009151_epn134/001/AO2D.root > log_o2.txt
# sometimes you may need this workflow: o2-analysis-collision-converter $CONFIG | 

# for running on file on alien, use alien:///alice/data/2022/LHC22t/529552/apass3/0800/o2_ctf_run00529552_orbit0007634688_tf0000009151_epn134/001/AO2D.root

# check the exit code, if it's non-zero then let the user know.
if [ $? -ne 0 ]; then
    echo "The O2 workflow exited with non-zero exit code, something probably went wrong! Please check the log: log_o2.txt."
fi

# let the user know the script is finished
notify-send "Your O2 analysis is done!" # pop-up
echo -e "\a" # sound notification
