#!/bin/zsh
# Step 3 codon-column tables (codon_columns.py from recheck4, unchanged), run sequentially (the script rewrites <T>.tagsense_orfs.faa,
# so never run two with the same T at once). Usage: zsh run_cols.sh <dir> threads
R=$1; TH=$2; P="nice -n 10 uv run --no-project --quiet python -I $R/scripts/codon_columns.py $R"
T=GCA_056209055.1; C1=GCA_934647165.1; C2=GCA_030529995.1; C3=GCA_963656345.1; J=GCA_030526225.1; S=GCA_045599485.1
for p in "$T $C1" "$T $C2" "$T $C3" "$T $S" "$T $J" "$J $C1" "$J $C2" "$J $C3" "$C1 $C2" "$C1 $C3" "$C2 $C1" "$C2 $C3" "$C3 $C1" "$C3 $C2" "$S $T"; do
  ${=P} ${=p} $TH; echo "$(date)"
done
