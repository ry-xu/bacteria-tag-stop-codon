#!/bin/zsh
# Step 4: tRNAscan-SE 2 (-B, --detail, -s isotype scores), 3 genomes at a time x 1 thread. Usage: zsh run_trnascan_par.sh <dir>
R=$1; S=$R/scripts
zsh $S/run_trnascan.sh $R GCA_937907055.1 GCA_017577285.1 GCA_934647165.1 &
zsh $S/run_trnascan.sh $R GCA_045628485.1 GCA_937886445.1 GCA_030529995.1 &
zsh $S/run_trnascan.sh $R GCA_030526225.1 GCA_045617125.1 &
wait
