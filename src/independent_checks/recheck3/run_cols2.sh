#!/bin/zsh
# Second batch: redo A2->A1 (first attempt's ORF query file was overwritten mid-run), plus control->control calibration
# pairs (TAG-as-stop genomes processed exactly like the tests). Usage: zsh run_cols2.sh <dir>
R=$1; P="nice -n 10 uv run --no-project --quiet python -I $R/scripts/codon_columns.py $R"
rm -f $R/out/cols/GCA_045628485.1__GCA_937907055.1.tblastn.tsv
for p in "GCA_045628485.1 GCA_937907055.1" "GCA_045617125.1 GCA_017577285.1" "GCA_045617125.1 GCA_937886445.1" "GCA_030529995.1 GCA_934647165.1"; do
  ${=P} ${=p} 3; echo "$(date)"
done
