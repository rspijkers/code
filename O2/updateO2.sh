# check if we need to build only one executable, or all executables in a directory:
while getopts ":d:e:" option; do
    case $option in
        d) # build a dir
            EXECUTABLE="$OPTARG/all";;
        e) # build an executable
            EXECUTABLE="O2Physicsexe-analysis-$OPTARG";;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

EXECDIR=$(pwd)
LOGFILE=log_updateO2.txt

# Enter installation directory 
cd $ALIBUILD_WORK_DIR/BUILD/O2Physics-latest/O2Physics/
if [ $? -ne 0 ]; then
    echo "Changing directory to O2Physics build directory failed, exiting!"
    exit
fi

# Redirect stdout and stderr to log file
exec > "$EXECDIR/$LOGFILE"
exec 2>&1
# Build
alienv setenv O2Physics/latest ninja/latest -c ninja install $EXECUTABLE
# reset std out to terminal
exec &>/dev/tty 

# Automatically check for any warnings/errors in the log file
ERRORS=`cat $EXECDIR/$LOGFILE | egrep -i 'error|warn|fatal'`
if [ -z "$ERRORS" ]; then
  echo "Build finished without encountering any warnings or errors"
else
  echo "Build finished with the following warnings and/or errors:"
  echo "$ERRORS" # double quotes are needed to keep linebreaks
fi