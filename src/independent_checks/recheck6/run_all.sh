#!/bin/zsh
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}  # repository root (portable copy)
# Rerun everything in recheck6 (~15 min on a loaded 12-core mini; < 2 GB RAM). Outputs in out/.
set -e
cd "${R6:-$REPO/work/recheck6}"  # working dir holding scripts/ (a copy of this folder), data/, out/
P=(uv run --no-project --with pyrodigal==3.7.1 python -I)
G=(GCA_056209055.1 QTP_TA20190011_bin.237 QTP_TA20190013_bin.296 GCA_056207115.1 GCA_030529995.1 GCA_963656345.1 GCA_045599485.1)
CNGB=${CNGB_DIR:-$REPO/data/cngb_bins}/Ass_MAGs_QS50_n6684  # authors' MAG set, CNGB CNP0001390
mkdir -p data out/blastdb out/tbn out/rf out/trna
for g in $G; do gzip -dc $REPO/data/genomes/$g.fna.gz > data/$g.fna; makeblastdb -in data/$g.fna -dbtype nucl -out out/blastdb/$g -parse_seqids > /dev/null; done
# 1. same species? contamination?
for b in QTP_TA20190011_bin.237 QTP_TA20190013_bin.296 GCA_056207115.1; do
  minimap2 -cx asm20 -t 4 --secondary=no data/GCA_056209055.1.fna data/$b.fna 2>/dev/null > out/${b}_vs_candidate.paf
  $P scripts/ani_contigs.py data/$b.fna out/${b}_vs_candidate.paf ${b}_vs_cand
done
# 2. genes, tblastn of relatives' proteomes vs raw DNA, stops inside alignments; homology-free open-frame test
$P scripts/predict.py $G
./scripts/run_tblastn.sh
./scripts/run_stops.sh > out/stops/summary.txt
$P scripts/stops_in_hsps.py GCA_056209055.1_byQ1 GCA_056209055.1 GCA_030529995.1 50 > /dev/null
$P scripts/stops_in_hsps.py GCA_056209055.1_byQ2 GCA_056209055.1 GCA_963656345.1 50 > /dev/null
$P scripts/downstream_orf.py $G > out/downstream_orf.tsv
# 3. RF1 (DIAMOND vs Swiss-Prot RF set) and the bin.296 contig_000065 chimera (see README for the blastn commands)
diamond makedb --in $REPO/data/reference/rf.faa -d out/rf/rf --quiet
cat out/genes/*.faa > out/rf/all_prot.faa
diamond blastp -q out/rf/all_prot.faa -d out/rf/rf -o out/rf/rf_hits.tsv --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --max-target-seqs 50 --evalue 1e-20 --threads 4 --quiet --more-sensitive
./scripts/run_trnascan.sh
# 4. reads: coordinate check, polymorphic sites, cross-mapping from other species, TAGs in the new bins
cat $(ls $CNGB/*.fa | grep -v -E "TA20190005_bin.134.fa|TA20190011_bin.237.fa|TA20190013_bin.296.fa") | makeblastdb -dbtype nucl -out out/blastdb/cngb_other -title cngb_other -in - > /dev/null
$P scripts/crossmap_check.py > out/crossmap_summary.txt
$P scripts/read_sites_check.py > out/read_sites_check.txt
$P scripts/sites_in_bins.py > out/sites_in_bins.txt
# 5. meaning: Swiss-Prot columns at TAG sites vs random TGG codons
$P scripts/column_check.py GCA_056209055.1 cand_loose > out/columns/summary.txt
