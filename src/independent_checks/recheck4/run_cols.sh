#!/bin/zsh
# Step 3 codon-column tables (codon_columns.py: test ORFs with TAG->X, tblastn vs other genome's raw DNA), run sequentially
# (the script rewrites <T>.tagsense_orfs.faa, so never run two with the same T at once). Usage: zsh run_cols.sh <dir> threads
R=$1; TH=$2; P="nice -n 10 uv run --no-project --quiet python -I $R/scripts/codon_columns.py $R"
L1=GCA_056232285.1; L2=GCA_056261325.1; C1=GCA_045640465.1; C2=GCA_902803485.1; C3=GCA_017959345.1
for p in "$L1 $C1" "$L1 $C2" "$L1 $C3" "$L2 $C1" "$L2 $C2" "$L2 $C3" "$L1 $L2" "$L2 $L1" "$C1 $C2" "$C1 $C3" "$C2 $C1" "$C2 $C3"; do
  ${=P} ${=p} $TH; echo "$(date)"
done
