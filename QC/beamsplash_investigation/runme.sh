# run this while having the QualityControl env loaded
# this script will most likely output a bunch of 

# Set the info like run nr, pass, etc here.
YEAR="2023"
RUNNR="537903"
PASS="apass1"
PERIOD="LHC23zc"
INTERVAL="1000" # paste the result of the first loop here

# names of the input and output files, hardcoded!
INFILE="QCfiles.dat"
OUTFILE="problematic_ctfs.dat"
if [[ -f "$INFILE" ]]; then
    rm $INFILE
fi

# make a list of the QC files per 10-minute time interval
for path in `alien_find /alice/data/$YEAR/$PERIOD/$RUNNR/$PASS/*/QC/001/QC.root`
do
    echo alien://${path} >> $INFILE
done
# run over them to find the problematic time interval(s)
root readQC_loop.C -q

# print all the problematic time intervals, in case there are more than 1 you have to handle them manually I guess
while read -r path
do
    echo $path
done < $OUTFILE

# extract the time interval from the path, probably only works if there is exactly one. 
INTERVAL=$(cut -d/ -f10 $OUTFILE)
echo $INTERVAL

# delete and make filelist again, this time with QC files per CTF
if [[ -f "$INFILE" ]]; then
    rm $INFILE
fi
for path in `alien_find /alice/data/$YEAR/$PERIOD/$RUNNR/$PASS/$INTERVAL/o2*/QC.root`
do
    echo alien://${path} >> $INFILE
done
# again find the problematic files
root readQC_loop.C -q

# Finally we are ready to check the time!
while read -r QC
do  
    # find the corresponding ctf's from the QC.root files:
    CTF=${QC//"$PASS"/"raw"} # replace "apass1" with "raw" to get the raw datafiles (CTF's)
    CTF=${CTF//"/QC"/} # remove the "/QC" at the end, so we don't get [ctf]/QC.root but [ctf].root
    # echo $CTF
    # now we print the time at which the raw datafile was last changed
    echo $(alien_stat "$CTF" | grep "Last change")
done < $OUTFILE