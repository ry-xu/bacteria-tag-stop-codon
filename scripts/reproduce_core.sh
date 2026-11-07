#!/usr/bin/env bash
# Reproduce the core per-genome results for the 12 TAG-reading genomes (about 10-20 min on a laptop, plus downloads).
# Prereqs: scripts/setup_reference_data.sh (Swiss-Prot, GTDB files); scripts/fetch_cngb_bins.sh (the two CNGB bins);
# diamond on PATH. NCBI genomes are downloaded on first use into data/genomes/ (or set GENOME_DIR to an existing cache).
# Outputs go to results/ (gitignored). Compare with the shipped tables in data/.
set -euo pipefail
cd "$(dirname "$0")/.."
READERS=$(cut -f1 data/genomes.tsv | tail -n +2)
mkdir -p results/meaning
# 1. read-through test per genome -> results/rt/<acc>.tsv + .stats.json   (compare: data/readthrough_counts.tsv)
for a in $READERS; do uv run python src/readthrough.py "$a"; done
# 2. release factors and stop usage -> results/rf_profile.tsv              (compare: data/rf1/rf_profile.tsv)
cut -f1 data/genomes.tsv | tail -n +2 > results/readers.txt
uv run python src/rf_profile.py results/readers.txt 4
# 3. RF1 neighbourhood (genes called with TAG as sense) -> results/locus/  (compare: data/rf1/locus_neighbourhoods.tsv)
uv run python src/rf1_locus.py GCA_937907055.1 GCA_045628485.1 GCA_030526225.1 GCA_056232285.1 GCA_056261325.1 GCA_056209055.1
# 4. what TAG means: relatives' residue at each TAG site, and the calibrated test (compare: data/tag_sites/, data/meaning/)
tail -n +2 data/tag_sites/relatives.tsv | while IFS=$'\t' read -r target lineage rels minid note; do
  uv run python src/ortholog_codons.py "$target" TAG 15 ${rels//,/ } > "results/meaning/${target}_codons.txt"
  MINID=$minid uv run python src/ortholog_profiles.py "$target" TAG 15 ${rels//,/ } > "results/meaning/${target}_profiles.txt"
done
# 5. screen ranking and the overview figure, from the shipped screen table (no per-genome files needed)
rm -rf results/rt_readers; mv results/rt results/rt_readers   # so rank_readers/figure2 read data/screen/screened_genomes.tsv.gz, not just these 12
uv run python src/rank_readers.py > results/rank_readers.txt
uv run python src/figure2.py
echo "done; see results/"
