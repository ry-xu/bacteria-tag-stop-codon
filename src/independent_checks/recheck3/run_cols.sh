#!/bin/zsh
# Step 3 codon-column tables (test ORFs with TAG->X, tblastn vs other genome). Usage: zsh run_cols.sh <dir>
R=$1; P="nice -n 10 uv run --no-project --quiet python -I $R/scripts/codon_columns.py $R"
for p in "GCA_937907055.1 GCA_017577285.1" "GCA_937907055.1 GCA_937886445.1" "GCA_937907055.1 GCA_045617125.1" \
         "GCA_045628485.1 GCA_017577285.1" "GCA_045628485.1 GCA_937886445.1" "GCA_045628485.1 GCA_045617125.1" \
         "GCA_030526225.1 GCA_934647165.1" "GCA_030526225.1 GCA_030529995.1" \
         "GCA_937907055.1 GCA_045628485.1" "GCA_045628485.1 GCA_937907055.1" \
         "GCA_030526225.1 GCA_937907055.1" "GCA_030526225.1 GCA_045628485.1"; do
  ${=P} ${=p} 3; echo "$(date)"
done
