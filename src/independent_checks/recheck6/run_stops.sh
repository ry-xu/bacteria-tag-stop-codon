#!/bin/zsh
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}  # repository root (portable copy)
# Stops inside high-identity tblastn alignments, per genome; strict (HSP >= 50% id) and loose (>= 35%) filters.
cd "${R6:-$REPO/work/recheck6}"  # working dir holding scripts/ (a copy of this folder), data/, out/
R() { uv run --no-project python -I scripts/stops_in_hsps.py "$@" }
A=GCA_030529995.1; Bq=GCA_963656345.1; C=GCA_056207115.1
for minpid in 50 35; do
  s=$([[ $minpid == 50 ]] && echo all || echo loose)
  R cand_$s GCA_056209055.1 $A,$Bq,$C $minpid
  R b237_$s QTP_TA20190011_bin.237 $A,$Bq,$C $minpid
  R b296_$s QTP_TA20190013_bin.296 $A,$Bq,$C $minpid
  R glu_$s GCA_045599485.1 $A,$Bq,$C $minpid
  R relB_$s $C $A,$Bq $minpid
  R caku1_$s $A $Bq,$C $minpid
  R caku2_$s $Bq $A,$C $minpid
done
# same single query proteome (GCA_056207115.1, the closest relative) for every other genome
for t in GCA_056209055.1 QTP_TA20190011_bin.237 QTP_TA20190013_bin.296 GCA_030529995.1 GCA_963656345.1 GCA_045599485.1; do
  R ${t}_byRelB $t $C 50
done
# CAKUFA01 queries only (both), for every non-CAKUFA01 genome
for t in GCA_056209055.1 QTP_TA20190011_bin.237 QTP_TA20190013_bin.296 GCA_056207115.1 GCA_045599485.1; do
  R ${t}_byCAKU $t $A,$Bq 50
done
