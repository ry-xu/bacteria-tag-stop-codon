#!/usr/bin/env bash
# The two lineage-6 genomes QTP_TA20190011_bin.237 and QTP_TA20190013_bin.296 are not in NCBI. They come from the
# authors' full MAG set of the Tibetan-plateau catalogue (Li et al., Microbiome 2025, 13:236), CNGB project CNP0001390,
# file Ass_MAGs_QS50_n6684.tar.gz. This streams the archive once (resumable byte ranges; the CNGB server drops long
# transfers and sometimes refuses connections for hours) and keeps only the two bins, saved under the names the scripts
# use: data/genomes/QTP_TA20190011_bin.237.fna.gz and data/genomes/QTP_TA20190013_bin.296.fna.gz.
# Uses bsdtar (macOS default tar); with GNU tar replace --include '<pat>' by --wildcards '<pat>'.
set -euo pipefail
cd "$(dirname "$0")/.."
URL=${CNGB_URL:-"https://ftp.cngb.org/pub/CNSA/data4/CNP0001390/Ass_MAGs_QS50_n6684.tar.gz"}
mkdir -p data/genomes tmp_cngb
uv run python src/range_stream.py "$URL" 32 \
  | tar -xzf - -C tmp_cngb --include '*TA20190011_bin.237.fa' --include '*TA20190013_bin.296.fa'
gzip -c "$(find tmp_cngb -name 'T1900044956_TA20190011_bin.237.fa')" > data/genomes/QTP_TA20190011_bin.237.fna.gz
gzip -c "$(find tmp_cngb -name 'T1900044957_TA20190013_bin.296.fa')" > data/genomes/QTP_TA20190013_bin.296.fna.gz
rm -rf tmp_cngb
ls -la data/genomes/QTP_*
