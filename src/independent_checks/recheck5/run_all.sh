#!/bin/zsh
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}  # repository root (portable copy)
# recheck5 (Aristaeellaceae MAG GCA_056209055.1, Tibetan wild ass faeces; claim: UAG read as Trp) - full reproduction.
# Usage: zsh run_all.sh <recheck5_dir>. Needs BLAST+, DIAMOND, uv (pyhmmer, pyrodigal, biopython), tRNAscan-SE (on PATH or $TRNASCAN_BIN).
# Inputs (no downloads): recheck4/data/{queries.faa,rf_queries.faa,rf12_sprot.hmm,gat_aars.faa}; data/sprot.dmnd symlinked.
# Scripts copied from recheck4 unless marked NEW (recheck5). Total runtime ~35 min at 3 threads.
set -e
R=$1; S=$R/scripts; G=$REPO/data/genomes
T=GCA_056209055.1; C1=GCA_934647165.1; C2=GCA_030529995.1; C3=GCA_963656345.1; J=GCA_030526225.1; SF=GCA_045599485.1
EC=GCF_000005845.2; BS=GCF_000009045.1; ALL="$T $C1 $C2 $C3 $J $SF $EC $BS"; CT=$C1,$C2,$C3; OG=$J,$SF,$EC,$BS
PY="nice -n 10 uv run --no-project --quiet"
mkdir -p $R/data $R/out/{rf,trna,cols,genes,dup,nuc,rt,rna}
cp $REPO/src/independent_checks/inputs/{queries.faa,rf_queries.faa,rf12_sprot.hmm,gat_aars.faa} $R/data/; ln -sf $REPO/data/sprot.dmnd $R/data/sprot.dmnd
for a in ${=ALL}; do gunzip -c $G/$a.fna.gz > $R/data/$a.fna; makeblastdb -in $R/data/$a.fna -dbtype nucl -out $R/data/$a -parse_seqids >/dev/null; done
# ---- step 1: stops inside conserved-protein alignments, spread, contamination
zsh $S/run_tblastn.sh $R ${=ALL}
python3 -I $S/stops_in_alignments.py $R ${=ALL}                                     # -> out/summary_stops.tsv, out/<acc>.stop_sites.tsv
python3 -I $S/contig_spread.py $R $T $J $SF; python3 -I $S/contig_perm.py $R $T
python3 -I $S/contam_check.py $R $T $C1 $C2 $C3 $J $SF
python3 -I $S/marker_inventory.py $R $T $CT                                         # NEW: families missing / extra copies in T
python3 -I $S/dup_copies.py $R $T dnaK,hisS,rnc,rph,rpmJ,rpsT,pyk,purA,glmS,guaA,lepA,aspS,map,atpA $CT,$J,$SF   # NEW
python3 -I $S/gc_profile.py $R $T $C1 $C2 $C3 $J $SF -- JBAEHL010000065.1           # NEW
for O in $SF $C1 $C2 $C3 $J; do nice -n 10 blastn -task dc-megablast -query $R/data/$T.fna -db $R/data/$O -evalue 1e-20 -num_threads 1 -max_target_seqs 5 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen" -out $R/out/nuc/T__$O.blastn.tsv; done   # shared DNA with SFMI01
for a in $T $C1 $C2 $C3 $J $SF; do nice -n 10 blastn -task dc-megablast -query $R/data/$a.fna -db $R/data/$a -evalue 1e-20 -num_threads 1 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" -out $R/out/dup/$a.self.tsv; python3 -I $S/self_dup.py $R $a 80 97; done   # NEW
# ---- step 2: release factors
$PY --with pyhmmer --with biopython python -I $S/rf_hmm.py $R 30 ${=ALL}                 # -> out/rf_hmm.tsv
python3 -I $S/rf_fragments.py $R $T $J $C1 > $R/out/rf/rf_fragments_tests.txt          # recheck4-fixed copy (min e-value)
$PY --with pyrodigal python -I $S/prfA_synteny.py $R 8 $T $C1 $C2 $C3 $J $SF          # recheck5: targets on the command line
$PY --with pyrodigal python -I $S/prfB_synteny.py $R 8 $T $C1 $J $SF                   # NEW (RF2 locus)
$PY --with biopython python -I $S/rf_compare.py $R > $R/out/rf/rf_compare.txt          # NEW: Swiss-Prot reciprocal, identities, motifs
python3 -I $S/rf1_rate.py $R $T $C2 $C3                                                # NEW: RF1 divergence percentile
$PY --with pyrodigal python -I $S/locus_genes.py $R T_tRNA:$T:JBAEHL010000072.1:10000:17979 T_prfA:$T:JBAEHL010000065.1:29000:36493
# ---- step 3: what does TAG encode
zsh $S/run_cols.sh $R 3; zsh $S/run_cols2.sh $R 3                                      # codon-column tables (T, J, controls, outgroups)
$PY --with pyrodigal python -I $S/ctrl_genes.py $R $T $C1 $C2 $C3 $J $SF $EC $BS        # pyrodigal genes (coding filter, read-through)
python3 -I $S/tag_meaning5.py $R $T $CT $OG coding > $R/out/tag_meaning5_T_coding.txt  # NEW: main call (+ W-column conservation)
python3 -I $S/tag_meaning5.py $R $T $CT $OG strict > $R/out/tag_meaning5_T_strict.txt
python3 -I $S/tag_meaning5.py $R $J $CT '' coding                                      # positive control (Glu reader)
for t in "$C1 $C2,$C3" "$C2 $C1,$C3" "$C3 $C1,$C2"; do python3 -I $S/tag_meaning5.py $R ${=t} '' coding; done   # negative controls
python3 -I $S/tag_meaning5_min1.py $R $T $CT; python3 -I $S/codon_by_column.py $R $T $CT; python3 -I $S/codon_by_column.py $R $J $CT
python3 -I $S/column_conservation.py $R $T $CT $OG
python3 -I $S/tag_site_table.py $R $T $CT $OG > $R/out/tag_sites_T.tsv
python3 -I $S/contig_origin.py $R $T $C1 $C2,$C3,$SF,$J
for t in "$T $CT" "$C1 $C2,$C3" "$J $CT"; do $PY --with pyrodigal python -I $S/readthrough.py $R ${=t}; done   # NEW: read-through by stop type
# ---- step 4: tRNA
zsh $S/run_trnascan.sh $R $T $C1 $C2 $C3 $J $SF $EC $BS
$PY --with biopython python -I $S/trna_trp_cta.py $R $T $C1 $C2 $C3 $J $SF $EC $BS > $R/out/trna/trna_trp_cta.txt   # NEW
