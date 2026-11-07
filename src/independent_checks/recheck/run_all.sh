#!/bin/zsh
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}  # repository root (portable copy)
# Reproduce the independent stop-codon / RF1 recheck. Usage: zsh run_all.sh <recheck_dir>
# Tools: BLAST+ 2.17.0 (Homebrew), DIAMOND 2.2.8, pyhmmer + biopython via uv, system python3.
set -e
R=$1; S=$R/scripts; SRC=$REPO/data
mkdir -p $R/data $R/out $R/dl/hmm
for a in GCA_949290385.1 GCA_965332555.1 GCA_045589595.1 GCA_015061915.1 GCA_934644965.1 GCF_000005845.2 GCF_000009045.1; do
  cp $SRC/genomes/$a.fna.gz $R/data/ && gunzip -f $R/data/$a.fna.gz
  makeblastdb -in $R/data/$a.fna -dbtype nucl -out $R/data/$a -parse_seqids >/dev/null
done
# 1. queries: 456 E. coli K-12 + B. subtilis 168 conserved proteins (Swiss-Prot, by gene name) + RF1/RF2/ArfB
python3 -I $S/extract_queries.py $SRC/uniprot_sprot.fasta.gz $R/data/queries.faa
#    rf_queries.faa = P0A7I0 (Ec prfA), P07012 (Ec prfB), P40711 (Ec arfB), P45872 (Bs prfA), P28367 (Bs prfB), pulled from the same Swiss-Prot file
# 2. tblastn (standard code, evalue 1e-5, 3 threads) of all queries vs every genome
zsh $S/run_tblastn.sh $R
# 3. in-alignment stop census (e<1e-10, qcov>=50%, codons deduplicated by genomic position)
python3 -I $S/stops_in_alignments.py $R GCA_949290385.1 GCA_965332555.1 GCA_045589595.1 GCA_015061915.1 GCA_934644965.1 GCF_000005845.2 GCF_000009045.1
python3 -I $S/contig_spread.py $R GCA_949290385.1 GCA_965332555.1
python3 -I $S/tag_site_orthologs.py $R GCA_949290385.1 GCA_965332555.1 GCA_045589595.1 GCA_015061915.1 GCA_934644965.1
python3 -I $S/tag_site_orthologs.py $R GCA_965332555.1 GCA_949290385.1 GCA_045589595.1 GCA_015061915.1 GCA_934644965.1
python3 -I $S/terminal_tag_check.py $R GCA_949290385.1 GCA_965332555.1 GCA_934644965.1 GCA_015061915.1
# 4. RF1/RF2: HMMs from InterPro API (PF03462, PF00472 Pfam; TIGR00019, TIGR00020 NCBIfam), 6-frame ORFs in std and TAG-sense codes
for p in PF03462 PF00472; do curl -sSL -o $R/dl/hmm/$p.hmm.gz "https://www.ebi.ac.uk/interpro/api/entry/pfam/$p?annotation=hmm"; done
for t in TIGR00019 TIGR00020; do curl -sSL -o $R/dl/hmm/$t.hmm.gz "https://www.ebi.ac.uk/interpro/api/entry/ncbifam/$t?annotation=hmm"; done
gunzip -f $R/dl/hmm/*.gz
uv run --quiet --with pyhmmer --with biopython python -I $S/rf_hmm.py $R $R/dl/hmm GCA_949290385.1 GCA_965332555.1 GCA_045589595.1 GCA_015061915.1 GCA_934644965.1 GCF_000005845.2 GCF_000009045.1
diamond blastp --db $SRC/sprot.dmnd --query $R/out/rf_candidates.faa --out $R/out/rf_candidates.sprot.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle --max-target-seqs 5 --threads 3 --sensitive --quiet
# 5. prfA synteny: control prfA +-15 kb ORFs -> tblastn vs test genomes; blastx of candidate regions vs Swiss-Prot
python3 -I $S/prfA_neighborhood.py $R
for t in GCA_949290385.1 GCA_965332555.1; do
  tblastn -query $R/data/prfA_neighbors.faa -db $R/data/$t -evalue 1e-10 -max_target_seqs 3 -num_threads 3 \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen" > $R/out/prfA_neighbors_vs_$t.tsv
done
