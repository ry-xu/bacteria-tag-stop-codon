#!/bin/bash
# Second pass over ERR14999989 (kangaroo sample), streamed from SRA, N streams of PER spots each:
#  - tloop_<i>.fq : read pairs where either mate has a tRNA T-loop motif (GTTC[GA]AATC or revcomp), for ARAGORN
#  - mag_<i>.sam  : reads mapping to the kangaroo MAG (hits only), for TAG-site read support and mate linkage
REPO=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p "$REPO/results/trna_hunt" && cd "$REPO/results/trna_hunt"   # outputs here; read by src/hunt_anticodons.py / hunt2_analyze.py
cp -n "$REPO/data/trna_hunt/clade_trnas.fa" "$REPO/data/trna_hunt/clade_trnas.tsv" . 2>/dev/null || true
[ -s kang_mag.mmi ] || { gzip -dc "$REPO/data/genomes/GCA_965332555.1.fna.gz" > kang_mag.fa && minimap2 -x sr -d kang_mag.mmi kang_mag.fa; }
N=${N:-4}; PER=${PER:-31250000}
for i in $(seq 0 $((N-1))); do
  a=$((i*PER+1)); b=$(((i+1)*PER))
  (fastq-dump --stdout --split-spot --skip-technical -N $a -X $b ERR14999989 2>fqd2_$i.err \
     | tee >(minimap2 -x sr -a -t 1 --sam-hit-only kang_mag.mmi - > mag_$i.sam 2>mm2_$i.err) \
     | awk 'NR%8==1{h1=$0} NR%8==2{s1=$0} NR%8==4{q1=$0} NR%8==5{h2=$0} NR%8==6{s2=$0}
            NR%8==0{ if (s1 ~ /GTTC[GA]AATC|GATT[TC]GAAC/ || s2 ~ /GTTC[GA]AATC|GATT[TC]GAAC/) print h1"\n"s1"\n+\n"q1"\n"h2"\n"s2"\n+\n"$0 }' \
     > tloop_$i.fq; echo "part $i done $(date)") &
done
wait
echo all done $(date)
