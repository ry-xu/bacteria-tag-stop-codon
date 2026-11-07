#!/bin/zsh
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}  # repository root (portable copy)
cd "${R6:-$REPO/work/recheck6}"  # working dir holding scripts/ (a copy of this folder), data/, out/
for g in GCA_056209055.1 QTP_TA20190011_bin.237 QTP_TA20190013_bin.296 GCA_056207115.1; do
  [[ -s out/trna/$g.tsv ]] && continue
  ${TRNASCAN_BIN:+$TRNASCAN_BIN/}tRNAscan-SE -B -q --thread 2 -o out/trna/$g.tsv -f out/trna/$g.ss --detail data/$g.fna > /dev/null 2>&1
  echo "$g $(date +%T)"
done
echo TRNADONE
