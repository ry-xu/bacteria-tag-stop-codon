#!/bin/zsh
# Positive control: E. coli selenoproteins (Sec = TGA) vs E. coli genome; census should report in-alignment TGA.
R=$1
tblastn -query $R/data/sec_queries.faa -db $R/data/GCF_000005845.2 -evalue 1e-10 -max_target_seqs 5 -num_threads 3 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sframe qseq sseq" \
  | awk -F'\t' 'BEGIN{OFS="\t"} {$1="Ec|sec_"substr($1,4,6)"|x"; print}' > $R/out/sec_control.tblastn.tsv
