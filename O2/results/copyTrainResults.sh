# this shell script downloads the results from a hyperloop train run to stbc. 
# The script expects a train number, and tries to get the test dir using this number. 
# Once the train numbers exceed 6 digits this procedure will fail, but hopefully I'll have finished my PhD by then
# Run it as follows from inside an O2 environment:
#   `bash copyTrainResults.sh 123456`
# or in one go like:
#   `alienv setenv O2Physics/latest -c bash copyTrainResults.sh 123456`

TRAINRUN=$1
# try to get the testdir from the train run number
DIGITS=${TRAINRUN:0:2}
TESTID=00$DIGITS/00$TRAINRUN

# make sure we are in the right directory...
cd /user/rspijker/project/code/O2/results

# make and enter the relevant dir
[ ! -d "$TRAINRUN" ] && mkdir $TRAINRUN
cd $TRAINRUN

### 1 Download the full merged output AnalysisResults.root
# there is an extra nested directory
NEST=`alien_ls /alice/cern.ch/user/a/alihyperloop/outputs/$TRAINRUN`

alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/$TRAINRUN/$NEST/AnalysisResults.root file:./AnalisysResults.root

### 2 Download config, stdout from the test directory
curl -k --cert ~/.globus/usercert.pem --cert-type PEM --key ~/.globus/userkey.pem \
https://alimonitor.cern.ch/train-workdir/tests/$TESTID/stdout.log -o stdout.log \
https://alimonitor.cern.ch/train-workdir/tests/$TESTID/dpl-config.json -o dpl-config.json \
https://alimonitor.cern.ch/train-workdir/tests/$TESTID/configuration.json -o configuration.json

# check if curl was succesful. 
CURLEXIT=$?
[[ $CURLEXIT != 0 ]] && echo curl failed with exit code $CURLEXIT!!!

# check if the first line of the json file starts with a '}', if not then curl probably downloaded the text of a 404 error...
JSONCHECK=$(head -n 1 dpl-config.json)
[[ $JSONCHECK != "{" ]] && echo "WARNING!! The first line of 'dpl-config.json' is not '{', you've probably just downloaded the text of a 404 error..."
# this check can probably be nicer, using grep 404 or something...
