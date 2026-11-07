#!/usr/bin/env bash
# Download the public reference data the analysis needs into data/ (all gitignored; ~0.6 GB on disk).
#   - UniProtKB/Swiss-Prot (current release) -> data/sprot.dmnd   (the original run used the release current on 2026-11-04;
#     a newer release can shift individual hits slightly)
#   - GTDB R232 metadata and trees          -> data/{bac120,ar53}_metadata.tsv.gz, data/{bac120,ar53}_r232.tree.gz
#   - DIAMOND databases for the shipped marker sets in data/reference/ -> data/{rf,arfb,markers}.dmnd
# Genomes are fetched on demand from NCBI by src/readthrough.py (cached in data/genomes/). Needs: curl, diamond.
set -euo pipefail
cd "$(dirname "$0")/.."
D=data
GTDB=https://data.gtdb.ecogenomic.org/releases/release232/232.0
[ -s $D/uniprot_sprot.fasta.gz ] || curl -fL -o $D/uniprot_sprot.fasta.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
[ -s $D/sprot.dmnd ] || diamond makedb --in $D/uniprot_sprot.fasta.gz -d $D/sprot --quiet
for f in bac120_metadata ar53_metadata; do
  [ -s $D/$f.tsv.gz ] || curl -fL -o $D/$f.tsv.gz $GTDB/${f}_r232.tsv.gz
done
for f in bac120_r232.tree.gz ar53_r232.tree.gz; do
  [ -s $D/$f ] || curl -fL -o $D/$f $GTDB/$f
done
diamond makedb --in $D/reference/rf.faa -d $D/rf --quiet
diamond makedb --in $D/reference/arfb.faa -d $D/arfb --quiet
diamond makedb --in $D/reference/markers.faa.gz -d $D/markers --quiet
echo "reference data ready in $D/"
