#!/bin/zsh
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}  # repository root (portable copy)
# recheck3 (Aristaeellaceae ZC3RG03 / JAUNDF01 TAG claim) - full reproduction. Usage: zsh run_all.sh <recheck3_dir>
# Needs: BLAST+, DIAMOND (Homebrew), uv (pyhmmer, pyrodigal, biopython), tRNAscan-SE 2 (on PATH or $TRNASCAN_BIN). No downloads.
set -e
R=$1; S=$R/scripts; G=$REPO/data/genomes; SP=$REPO/data/uniprot_sprot.fasta.gz
T="GCA_937907055.1 GCA_045628485.1 GCA_030526225.1"; CA="GCA_017577285.1 GCA_937886445.1 GCA_045617125.1"; CB="GCA_934647165.1 GCA_030529995.1"; X="GCF_000005845.2 GCF_000009045.1"
PY="nice -n 10 uv run --no-project --quiet"
mkdir -p $R/data $R/out/rf $R/out/trna
for a in ${=T} ${=CA} ${=CB} ${=X}; do gunzip -c $G/$a.fna.gz > $R/data/$a.fna; makeblastdb -in $R/data/$a.fna -dbtype nucl -out $R/data/$a -parse_seqids >/dev/null; done
cp $REPO/src/independent_checks/inputs/{queries.faa,rf_queries.faa,gat_aars.faa,rf12_sprot.hmm} $R/data/  # portable copy: shipped outputs of extract_queries.py (needs recheck2 work files)
zsh $S/run_tblastn.sh $R ${=T} ${=CA} ${=CB} ${=X}                                     # step 1 tblastn
$PY python -I $S/stops_in_alignments.py $R ${=T} ${=CA} ${=CB} ${=X}                    # step 1 census -> out/summary_stops.tsv
$PY python -I $S/contig_spread.py $R ${=T}; $PY python -I $S/contig_perm.py $R ${=T}
$PY python -I $S/terminal_tag_check.py $R ${=T} ${=CA} ${=CB}; $PY python -I $S/frameshift_check.py $R ${=T} ${=CA} ${=CB} ${=X}
$PY python -I $S/contam_check.py $R ${=T} ${=CA} ${=CB} ${=X}
$PY python -I $S/contig_qc.py $R GCA_937907055.1 GCA_045628485.1,GCA_017577285.1,GCA_937886445.1,GCA_045617125.1
$PY python -I $S/contig_qc.py $R GCA_045628485.1 GCA_937907055.1,GCA_017577285.1,GCA_937886445.1,GCA_045617125.1
$PY python -I $S/contig_qc.py $R GCA_030526225.1 GCA_934647165.1,GCA_030529995.1
# step 2
$PY --with pyhmmer python -I $S/build_rf_hmms.py $R $SP                                 # RF1/RF2 HMMs from Swiss-Prot
diamond makedb --in $SP -d $R/data/sprot --threads 3 --quiet
$PY --with pyhmmer --with biopython python -I $S/rf_hmm.py $R 30 ${=T} ${=CA} ${=CB} ${=X}   # -> out/rf_hmm.tsv
$PY python -I $S/rf_fragments.py $R ${=T} GCA_934647165.1                              # RF1 fragments in raw DNA
$PY --with pyrodigal python -I $S/prfA_synteny.py $R 10                                 # -> out/rf/prfA_synteny.tsv
$PY --with pyrodigal python -I $S/locus_genes.py $R A1:GCA_937907055.1:CALCTW010000132.1:70000:83557 A2:GCA_045628485.1:JAZQPK010000024.1:1:8489 \
   cA_17577:GCA_017577285.1:PITV01000001.1:44000:58000 B_prmC:GCA_030526225.1:JAUNDF010000010.1:17000:30000 cB_29995:GCA_030529995.1:JAUMVT010000031.1:1:16000
# step 3
zsh $S/run_cols.sh $R; zsh $S/run_cols2.sh $R
$PY --with pyrodigal python -I $S/ctrl_genes.py $R ${=CA} ${=CB}
$PY python -I $S/tag_meaning.py $R GCA_937907055.1 GCA_017577285.1,GCA_937886445.1,GCA_045617125.1 > $R/out/tag_meaning_A1.txt
$PY python -I $S/tag_meaning.py $R GCA_045628485.1 GCA_017577285.1,GCA_937886445.1,GCA_045617125.1 > $R/out/tag_meaning_A2.txt
$PY python -I $S/tag_meaning.py $R GCA_030526225.1 GCA_934647165.1,GCA_030529995.1 > $R/out/tag_meaning_B.txt
$PY python -I $S/tag_meaning.py $R GCA_045617125.1 GCA_017577285.1,GCA_937886445.1 > $R/out/tag_meaning_ctrlA.txt   # calibration: TAG-stop genome
$PY python -I $S/tag_meaning.py $R GCA_030529995.1 GCA_934647165.1 > $R/out/tag_meaning_ctrlB.txt
$PY python -I $S/tag_swap.py $R GCA_937907055.1 GCA_045628485.1 GCA_017577285.1,GCA_937886445.1,GCA_045617125.1 > $R/out/tag_swap_A1_A2.txt
$PY python -I $S/tag_swap.py $R GCA_045628485.1 GCA_937907055.1 GCA_017577285.1,GCA_937886445.1,GCA_045617125.1 > $R/out/tag_swap_A2_A1.txt
$PY python -I $S/tag_swap.py $R GCA_030526225.1 GCA_937907055.1 GCA_934647165.1,GCA_030529995.1 > $R/out/tag_swap_B_GCA_937907055.1.txt
# step 4
zsh $S/run_trnascan_par.sh $R
$PY --with biopython python -I $S/trna_report.py $R ${=T} ${=CA} ${=CB} > $R/out/trna/trna_report.txt
$PY --with pyrodigal python -I $S/locus_genes.py $R trnaA1:GCA_937907055.1:CALCTW010000059.1:1500:14000 trnaA2:GCA_045628485.1:JAZQPK010000196.1:1:12500
