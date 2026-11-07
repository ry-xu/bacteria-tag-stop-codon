#!/bin/bash
# Stream kangaroo-sample reads (ERR14999989) from SRA, keep only reads mapping to the TAG-reading clade's tRNAs.
# 4 streams over spots 1..125M (~1/4 of the run); output: part_<i>.sam (hits only)
REPO=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p "$REPO/results/trna_hunt" && cd "$REPO/results/trna_hunt"   # outputs here; read by src/hunt_anticodons.py / hunt2_analyze.py
cp -n "$REPO/data/trna_hunt/clade_trnas.fa" "$REPO/data/trna_hunt/clade_trnas.tsv" . 2>/dev/null || true
for i in 0 1 2 3; do
  a=$((i*31250000+1)); b=$(((i+1)*31250000))
  (fastq-dump --stdout --split-spot --skip-technical -N $a -X $b ERR14999989 2>fqd_$i.err \
     | minimap2 -x sr -a -t 1 --sam-hit-only clade_trnas.fa - > part_$i.sam 2>mm_$i.err; echo "part $i done $(date)") &
done
wait
echo all done $(date)
