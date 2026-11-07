#!/bin/zsh
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}  # repository root (portable copy)
# recheck4 (Lachnospiraceae Tibetan-horse MAGs GCA_056232285.1 / GCA_056261325.1, TAG->Gln claim) - full reproduction.
# Usage: zsh run_all.sh <recheck4_dir>. Needs BLAST+, DIAMOND, uv (pyhmmer, pyrodigal, biopython), tRNAscan-SE (on PATH or $TRNASCAN_BIN).
# Inputs copied (no downloads): recheck3/data/{queries.faa,rf_queries.faa,rf12_sprot.hmm,gat_aars.faa}; data/sprot.dmnd symlinked.
# Scripts reused from recheck3 (some with thread counts lowered); rf_fragments.py has a recheck4 bug fix (locus filter max->min e-value).
set -e
R=$1; S=$R/scripts; G=$REPO/data/genomes
L1=GCA_056232285.1; L2=GCA_056261325.1; C1=GCA_045640465.1; C2=GCA_902803485.1; C3=GCA_017959345.1; X="GCF_000005845.2 GCF_000009045.1"
ALL="$L1 $L2 $C1 $C2 $C3 $X"; CT=$C1,$C2,$C3
PY="nice -n 10 uv run --no-project --quiet"
mkdir -p $R/data $R/out/{rf,trna,cols,genes,dup}
cp $REPO/src/independent_checks/inputs/{queries.faa,rf_queries.faa,rf12_sprot.hmm,gat_aars.faa} $R/data/; ln -sf $REPO/data/sprot.dmnd $R/data/sprot.dmnd
for a in ${=ALL}; do gunzip -c $G/$a.fna.gz > $R/data/$a.fna; makeblastdb -in $R/data/$a.fna -dbtype nucl -out $R/data/$a -parse_seqids >/dev/null; done
# step 1
zsh $S/run_tblastn.sh $R ${=ALL}
$PY python -I $S/stops_in_alignments.py $R ${=ALL}                  # -> out/summary_stops.tsv
$PY python -I $S/contig_spread.py $R $L1 $L2; $PY python -I $S/contig_perm.py $R $L1 $L2
$PY python -I $S/terminal_tag_check.py $R $L1 $L2 $C1 $C2 $C3; $PY python -I $S/frameshift_check.py $R ${=ALL}; $PY python -I $S/contam_check.py $R ${=ALL}
# step 2
$PY --with pyhmmer --with biopython python -I $S/rf_hmm.py $R 30 ${=ALL}        # -> out/rf_hmm.tsv, out/rf_candidates.faa
python3 -I $S/rf_fragments.py $R $L1 $L2 > $R/out/rf/rf_fragments_tests.txt     # control RF1/RF2 + Ec/Bs as queries
$PY --with pyrodigal python -I $S/prfA_synteny.py $R 10                          # -> out/rf/prfA_synteny.tsv
$PY --with pyrodigal python -I $S/locus_genes.py $R L1_prmC:$L1:JBADEI010000336.1:1:6094 L1_yqhQ:$L1:JBADEI010000307.1:1:3598 \
   L2_prmC:$L2:JBACEE010000374.1:1:3795 L2_yqhQ:$L2:JBACEE010000307.1:1:5058 C1_prfA:$C1:JAZUFC010000276.1:1:5000 C2_prfA:$C2:CADCQD010000057.1:1:7000 \
   L1_tRNA:$L1:JBADEI010000329.1:4000:13000 L2_c117:$L2:JBACEE010000117.1:1:3000
# step 3
zsh $S/run_cols.sh $R 2
$PY --with pyrodigal python -I $S/ctrl_genes.py $R $C1 $C2 $C3
$PY python -I $S/tag_meaning.py $R $L1 $CT > $R/out/tag_meaning_L1.txt; $PY python -I $S/tag_meaning.py $R $L2 $CT > $R/out/tag_meaning_L2.txt
$PY python -I $S/tag_meaning.py $R $C1 $C2,$C3 > $R/out/tag_meaning_ctrlC1.txt; $PY python -I $S/tag_meaning.py $R $C2 $C1,$C3 > $R/out/tag_meaning_ctrlC2.txt
python3 -I $S/tag_meaning4.py $R $L1 $L2 $CT > $R/out/tag_meaning4_L1.txt; python3 -I $S/tag_meaning4.py $R $L2 $L1 $CT > $R/out/tag_meaning4_L2.txt
python3 -I $S/tag_meaning4.py $R $C1 $C2 $C2,$C3 > $R/out/tag_meaning4_ctrlC1.txt; python3 -I $S/tag_meaning4.py $R $C2 $C1 $C1,$C3 > $R/out/tag_meaning4_ctrlC2.txt
python3 -I $S/tag_at_E.py $R $L2 $L1 $CT E; python3 -I $S/tag_at_E.py $R $L1 $L2 $CT E
# chimera / contamination
python3 -I $S/contig_origin.py $R $L2 $L1 $CT; python3 -I $S/contig_origin.py $R $L1 $L2 $CT
python3 -I $S/rp_copies.py $R $L1 $L2; python3 -I $S/paralog_ident.py $R 'Bs|rplK|NP_387983.2' $L2 $L1 $C1 $C2 $C3
# step 4/5
zsh $S/run_trnascan.sh $R $L1 $L2 $C1 $C2 $C3
$PY --with biopython python -I $S/trna_report.py $R $L1 $L2 $C1 $C2 $C3 > $R/out/trna/trna_report.txt
$PY --with biopython python -I $S/trna_compare.py $R $L1 $L2 $C1 $C2 $C3 > $R/out/trna/trna_compare.txt
