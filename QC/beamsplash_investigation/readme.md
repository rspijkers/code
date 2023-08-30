This folder has the tools necessary to correlate anomalies in phi-distributions with the real-time at which they occur. 
As these anomalies are usually caused by beam splash events, it is important to know the precise time of the anomaly. 
We will first look at the QC objects per 10-minute interval to find the approximate time. Then we will look per CTF-QC objects, and check the date of the corresponding raw ctf to see the exact time of the anomaly.

1. Put the relevant info in `runme.sh` (things like run number, pass, period, etc.)
2. Put a point in the problematic phi region in `readQC_loop.C`, such that we can use it as a proxy for problematic files. 
3. Run the script with `bash runme.sh`, it will output a bunch of info on accessing the files on the grid but the last lines are the timestamps of the problematic files. 
4. You can now correlate these timestamps with possible events like beam splashes! You can find info on these things in the logbook or in JIRA tickets. 