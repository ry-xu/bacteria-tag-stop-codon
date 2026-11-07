#!/bin/zsh
# Usage: run_tblastn.sh <recheck_dir>
R=$1
cat $R/data/queries.faa $R/data/rf_queries.faa > $R/data/all_queries.faa
for a in GCA_949290385.1 GCA_965332555.1 GCA_045589595.1 GCA_015061915.1 GCA_934644965.1 GCF_000005845.2 GCF_000009045.1; do
  nice -n 10 tblastn -query $R/data/all_queries.faa -db $R/data/$a -db_gencode 1 \
    -evalue 1e-5 -max_target_seqs 200 -num_threads 3 \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sframe qseq sseq" \
    -out $R/out/$a.tblastn.tsv
  echo "done $a $(date)"
done
