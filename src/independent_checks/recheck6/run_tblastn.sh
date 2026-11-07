#!/bin/zsh
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}  # repository root (portable copy)
# tblastn every query proteome (pyrodigal, table 11) against every other genome's raw DNA (code 11, stops shown as *)
cd "${R6:-$REPO/work/recheck6}"  # working dir holding scripts/ (a copy of this folder), data/, out/
Q=(GCA_030529995.1 GCA_963656345.1 GCA_056207115.1)
T=(GCA_056209055.1 QTP_TA20190011_bin.237 QTP_TA20190013_bin.296 GCA_056207115.1 GCA_030529995.1 GCA_963656345.1 GCA_045599485.1)
for q in $Q; do for t in $T; do
  [[ $q == $t ]] && continue
  o=out/tbn/${q}__${t}.tsv
  [[ -s $o ]] && continue
  tblastn -query out/genes/$q.faa -db out/blastdb/$t -db_gencode 11 -seg no -evalue 1e-10 -num_threads 4 -max_target_seqs 20 \
    -outfmt "6 qseqid saccver pident length qstart qend qlen sstart send sframe evalue bitscore qseq sseq" > $o.tmp && mv $o.tmp $o
  echo "$q $t $(wc -l < $o) $(date +%T)"
done; done
echo TBNDONE
