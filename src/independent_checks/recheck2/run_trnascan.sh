#!/bin/zsh
# Usage: zsh run_trnascan.sh <recheck_dir> acc...   (tRNAscan-SE 2.x, bacterial mode, isotype-specific models on)
R=$1; shift
export PATH=${TRNASCAN_BIN:+$TRNASCAN_BIN:}$PATH  # tRNAscan-SE 2: on PATH, or set TRNASCAN_BIN to its bin dir
for a in "$@"; do
  nice -n 10 tRNAscan-SE -B --thread 2 -q --detail -o $R/out/trna/$a.trna.txt -f $R/out/trna/$a.struct.txt \
    -s $R/out/trna/$a.iso.txt -a $R/out/trna/$a.trna.fa $R/data/$a.fna > $R/out/trna/$a.log 2>&1
  echo "done $a $(date)"
done
