# This script is meant to be executed from inside the O2 environment, perhaps this should be changed
# It aims to execute an O2 workflow, redirecting the std::out to log_o2.txt

# TODO: add option to pass txtfile with AO2D's to run over? maybe pass flags like '--mc' or '--run2' to automatically select the correct filelist?
# we can let the script check if we are on nikhef machines or not

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
HOST=`hostname`
if [ $HOST == stbc* ]; then # it seems we are on one of the stbc nodes at nikhef
    echo "we are on stbc"
elif [ $HOST == fiora ]; then
    echo "you are on the Nikhef login server, don't run O2 here!!"
    exit
fi
# if we are not on stbc or fiora, assume we are running locally
FILELIST="@Run3Samplefiles.txt"

# define the config so we can conveniently pass it to each workflow
CONFIG="-b --configuration json://ssbarconfig.json" # -b means no debug gui

# execute the entire workflow
o2-analysis-timestamp $CONFIG | o2-analysis-track-propagation $CONFIG | o2-analysis-event-selection $CONFIG | o2-analysis-multiplicity-table $CONFIG | o2-analysis-pid-tpc $CONFIG | o2-analysis-lf-lambdakzerobuilder $CONFIG | o2-analysis-lf-cascadebuilder $CONFIG | o2-analysis-lf-strangecorrelations $CONFIG --aod-file @Run3Samplefiles.txt > log_o2.txt

# check the exit code, if it's non-zero then let the user know.
if [ $? -ne 0 ]; then
    echo "The O2 workflow exited with non-zero exit code, something probably went wrong! Please check the log: log_o2.txt."
fi

# o2-analysis-timestamp $CONFIG | o2-analysis-track-propagation $CONFIG | o2-analysis-event-selection $CONFIG | o2-analysis-multiplicity-table $CONFIG | o2-analysis-pid-tpc $CONFIG | o2-analysis-lf-singlestrangebuilder $CONFIG | o2-analysis-lf-multistrangebuilder $CONFIG | o2-analysis-lf-strangecorrelations $CONFIG > log_o2.txt
