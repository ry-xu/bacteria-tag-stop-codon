#!/bin/zsh
# Outgroup codon-column tables: T vs E. coli / B. subtilis (for W-column conservation in tag_meaning5.py). Usage: zsh run_cols2.sh <dir> threads
R=$1; TH=$2; P="nice -n 10 uv run --no-project --quiet python -I $R/scripts/codon_columns.py $R"
T=GCA_056209055.1
for p in "$T GCF_000005845.2" "$T GCF_000009045.1"; do ${=P} ${=p} $TH; echo "$(date)"; done
