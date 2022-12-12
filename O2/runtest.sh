# Tutorial:
# o2-analysistutorial-histograms --aod-file $HOME/alice/data/AO2D-LHC21k6-1.root | o2-analysis-track-propagation | o2-analysis-timestamp > log_o2.txt


# lambdakzero:
#CONFIG="--configuration json://test.json"

# o2-analysis-timestamp $CONFIG | o2-analysis-track-propagation $CONFIG | o2-analysis-event-selection $CONFIG | o2-analysis-multiplicity-table $CONFIG | o2-analysis-pid-tpc $CONFIG | o2-analysis-lf-lambdakzerobuilder $CONFIG | o2-analysis-lf-lambdakzeroanalysis-mc --aod-file @aodfiles.txt $CONFIG


# now let's try cascade... seems it also works!
# specify aodfiles in json:
# CONFIG="-b --configuration json://cascadeconfig.json" # -b means no debug gui
# TODO: in cascade-builder, "d_UseWeightedPCA": "false" because of bug. setting should be true. 

# o2-analysis-timestamp $CONFIG | o2-analysis-track-propagation $CONFIG | o2-analysis-event-selection $CONFIG | o2-analysis-multiplicity-table $CONFIG | o2-analysis-pid-tpc $CONFIG | o2-analysis-lf-lambdakzerobuilder $CONFIG | o2-analysis-lf-lambdakzeroanalysis-mc $CONFIG | o2-analysis-lf-cascadebuilder $CONFIG | o2-analysis-lf-cascadeanalysismc $CONFIG > log_o2.txt


# Let's try now our custom analysis
CONFIG="-b --configuration json://ssbarconfig.json" # -b means no debug gui

o2-analysis-timestamp $CONFIG | o2-analysis-track-propagation $CONFIG | o2-analysis-event-selection $CONFIG | o2-analysis-multiplicity-table $CONFIG | o2-analysis-pid-tpc $CONFIG | o2-analysis-lf-lambdakzerobuilder $CONFIG | o2-analysis-lf-cascadebuilder $CONFIG | o2-analysis-lf-strangecorrelations $CONFIG > log_o2.txt
