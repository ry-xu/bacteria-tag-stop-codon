#!/bin/zsh
# Reproduce recheck2 (SFMI01 TAG reassignment). Usage: zsh run_all.sh <recheck2_dir>
# Tools: BLAST+ (Homebrew), pyhmmer/biopython/pyrodigal via uv, tRNAscan-SE 2 (on PATH or $TRNASCAN_BIN), python3.
set -e
R=$1; S=$R/scripts
T="GCA_024404275.1 GCA_045599485.1"; C="GCA_934796175.1 GCA_004556155.1 GCA_963624325.1"; X="GCF_000005845.2 GCF_000009045.1"
mkdir -p $R/dl $R/data $R/out
for a in ${=T} ${=C}; do mkdir -p $R/dl/$a; curl -sSf -o $R/dl/$a/pkg.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$a/download?include_annotation_type=GENOME_FASTA"; done
for a in ${=X}; do mkdir -p $R/dl/$a; curl -sSf -o $R/dl/$a/pkg.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$a/download?include_annotation_type=GENOME_FASTA,PROT_FASTA,GENOME_GFF"; done
for a in ${=T} ${=C} ${=X}; do (cd $R/dl/$a && unzip -o -q pkg.zip); cp $R/dl/$a/ncbi_dataset/data/$a/*_genomic.fna $R/data/$a.fna; makeblastdb -in $R/data/$a.fna -dbtype nucl -out $R/data/$a -parse_seqids >/dev/null; done
python3 -I $S/extract_queries.py $R                          # 457 Ec/Bs conserved proteins + RF1/RF2/ArfB/PrmC refs
zsh $S/run_tblastn.sh $R ${=T} ${=C} ${=X}                    # step 1 tblastn
python3 -I $S/stops_in_alignments.py $R ${=T} ${=C} ${=X}     # step 1 census (e<1e-10, qcov>=0.5, codon dedup)
python3 -I $S/contig_spread.py $R ${=T}; python3 -I $S/contig_perm.py $R ${=T}
python3 -I $S/terminal_tag_check.py $R ${=T} ${=C}
python3 -I $S/frameshift_check.py $R ${=T} ${=C} ${=X}; python3 -I $S/contam_check.py $R ${=T} ${=C} ${=X}
mkdir -p $R/out/ani; for q in ${=T} GCA_934796175.1; do for s in ${=T} ${=C}; do if [[ $q != $s ]]; then python3 -I $S/ani.py $R $q $s; fi; done; done
mkdir -p $R/dl/hmm   # HMMs: InterPro API pfam PF03462/PF00472, ncbifam TIGR00019/TIGR00020 (?annotation=hmm), gunzipped
uv run --quiet --with pyhmmer --with biopython python -I $S/rf_hmm.py $R ${=T} ${=C} ${=X}   # step 2
python3 -I $S/prfA_neighborhood.py $R; uv run --quiet --with pyrodigal python -I $S/prfA_synteny.py $R
for p in "GCA_024404275.1 GCA_934796175.1" "GCA_024404275.1 GCA_004556155.1" "GCA_024404275.1 GCA_963624325.1" "GCA_024404275.1 GCA_045599485.1" \
         "GCA_045599485.1 GCA_934796175.1" "GCA_045599485.1 GCA_004556155.1" "GCA_045599485.1 GCA_963624325.1" "GCA_045599485.1 GCA_024404275.1" \
         "GCA_963624325.1 GCA_934796175.1" "GCA_934796175.1 GCA_963624325.1"; do python3 -I $S/codon_columns.py $R ${=p} 3; done   # step 3
python3 -I $S/tag_meaning.py $R > $R/out/tag_meaning.txt
zsh $S/run_trnascan.sh $R ${=T} ${=C}                         # step 4
python3 -I $S/trna_extract.py $R ${=T} ${=C}
uv run --quiet --with biopython python -I $S/trna_compare.py $R
