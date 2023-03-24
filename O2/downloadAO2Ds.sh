# usage: enter O2 env, then do `bash downloadAO2D.sh`

# OUTPATH="${HOME}/alice/data/LHC22t_apass3/529552"
# [ ! -d "${OUTPATH}/$i/" ] && mkdir -p "${OUTPATH}/$i/"

# # this is for AO2D's that are not yet merged, so still in the o2_ctf_run...epn... format
# for i in {001..003}
# do
#     alien_cp alien:/alice/data/2022/LHC22t/529552/apass3/0800/o2_ctf_run00529552_orbit0007634688_tf0000009151_epn134/${i}/AO2D.root file:${OUTPATH}/${i}/AO2D.root 
# done

# Let's try to download a sample of LHC22t_triggersel to stbc
NFILES=4
RUNNR=529664
PERIOD="LHC22t"
PASS="apass3"

ALIPATH="/alice/data/2022/${PERIOD}/${RUNNR}/${PASS}"
OUTPATH="${HOME}/alice/data/${PERIOD}_${RUNNR}_${PASS}"

c=1
echo downloading the following files: 
for path in `alien_find -r ${ALIPATH} [0-9]{4}/o2_ctf.*/[0-9]{3}/AO2D\.root | head -${NFILES}`
do
    echo alien:${path} $c
    # alien_cp alien:${path} file:${OUTPATH}/AO2D-${c}.root
    ((c++))
done

