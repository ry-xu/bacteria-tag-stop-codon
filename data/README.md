# Data

Small tables only. Genomes are not included: NCBI accessions download on demand (`src/readthrough.py`), and the two
CNGB bins come from `scripts/fetch_cngb_bins.sh`. Reference databases come from `scripts/setup_reference_data.sh`.
Everything here was produced by the code in `src/` (or by the independent checks in `src/independent_checks/`). Each
section says which script made it.

`L1` to `L6` are the six lineages: L1 Victivallales, L2 SFMI01, L3 ZC3RG03, L4 JAUNDF01, L5 Lachnospiraceae near
UBA1267, L6 HGM11604. L2 to L4 and L6 are all in family Aristaeellaceae.

## The 12 genomes
- `genomes.tsv`: one row per TAG-reading genome. Columns: accession, lineage, taxonomy, source, host, place, date,
  BioProject/BioSample, completeness/contamination (with where the estimate comes from), contigs, size,
  read-through counts (`read-through genes / genes tested`), RF1 status, UAG-reading tRNA, inferred meaning, and two
  confidence columns. The first says how sure we are that TAG is read as sense; the second, how sure we are of the
  amino acid.
- `readthrough_counts.tsv`: per-genome counts from `src/readthrough.py` and `src/rf_profile.py`, for the 12 readers and
  the RF1-keeping relatives used as controls. `sampling_rates` is `{}` for full runs. Otherwise it gives the fraction
  of genes tested for each stop codon (screen runs).
- `readthrough_genes/<acc>.tsv`: every read-through gene of each reader, one row per gene. It gives the stop codon,
  the best Swiss-Prot subject, the subject position (`spos`) and residue (`sres`) aligned to the stop, and the
  alignment window around it (`qwin`/`swin`; `X` = the stop).

## What TAG means
- `tag_sites/<acc>_TAG.tsv`: one row per in-frame TAG site x relative. Columns: gene (`contig|begin|end|strand`,
  called with TAG as sense, table 15), amino-acid position, relative, percent identity of the protein alignment, and
  the relative's codon and residue at that position. Made by `src/ortholog_codons.py` (identity >= 50%). Relatives
  are in `tag_sites/relatives.tsv`. `_TAG.summary.txt` is the script's printed summary. All 12 were regenerated for
  this repository. Where an original output exists (JAUNDF01), the rerun is identical.
- `meaning/`: outputs of the calibrated test `src/ortholog_profiles.py` and the conservation-strata test
  `src/codon_strata.py`. `meaning/INDEX.tsv` lists each file with its controls and identity cut-off.

## RF1 (release factor 1, prfA)
- `rf1/breakpoints.tsv`: the deletion junctions, with coordinates and the gene order in an RF1-keeping relative.
- `rf1/locus_neighbourhoods.tsv`: genes around bL31, prmC, lepA, topIII and RF1 in readers and relatives, each with
  its best Swiss-Prot hit (`src/rf1_locus.py`).
- `rf1/rf_profile.tsv`: RF1/RF2 hits and stop-codon usage for 953 genomes (`src/rf_profile.py`). `RF1_motif` /
  `RF2_motif` = number of hits carrying the PxT / SPF stop-recognition motif.

## tRNAs
- `trnas/cua_trnas.fa`: the CUA-anticodon (UAG-reading) tRNAs, with coordinates in each genome. L6's tRNA-Trp is
  included because the L6 CUA tRNA is a copy of it.
- `trnas/cua_trna_isotype_scores.tsv`: tRNAscan-SE 2 isotype-specific model scores for the L2, L3 and L5 tRNAs.
- `trnas/trnascan_GCA_056209055.1.*`: tRNAscan-SE output for L6. The CUA tRNA is called "Undet".
- `trna_hunt/`: the L1 search for a UAG-reading tRNA in raw kangaroo-sample reads (ERR14999989). Included: the
  clade's own tRNAs, the CTA-anticodon tRNAs found in reads, and summaries.

## Raw-read support
- `read_support/`: per-run summaries (`*.txt`) and per-site codon tallies (`*.sites.tsv`: contig, 0-based position,
  strand, codons seen in reads) from `src/read_support.py`. `pooled_*.txt` pools L6 across four animals
  (`src/read_support_pool.py`). `_comp` = reads mapped competitively to the L6 genome and the L2 genome from the same
  animal. Runs: SRR17543912 (L2 horse), SRR6425796 (L3 horse), SRR24915835 (L4 warthog), CNCB GSA CRA013436 runs
  CRR939056/58/60/61 (L6 wild asses). L1 kangaroo read support is in `trna_hunt/hunt2_summary.txt`.

## Independent re-checks
- `independent_checks/recheck*/`: key outputs of six separate checks. Each was run by a different agent that did not
  read the main pipeline's code or results. Each used its own method: tblastn of ~456 conserved E. coli / B. subtilis
  proteins, or of relatives' proteomes, against raw genome DNA, with no gene calls. Coverage: recheck = L1, recheck2 =
  L2, recheck3 = L3 + L4, recheck4 = L5, recheck5 + recheck6 = L6. `summary_stops*.tsv` counts stops inside
  alignments, and `tag_meaning*.txt` gives the residue at TAG sites in the controls. Their scripts are in
  `src/independent_checks/`.

## The screen
- `screen/screened_genomes.tsv.gz`: all 42,230 genome accessions tested: per-stop genes tested, read-through counts,
  sampling rates, and the screen set(s) each genome came from. About 41,200 bacterial genomes after merging
  GenBank/RefSeq copies; the rest are ~900 archaea used as controls. 20 accessions have no set label: controls and
  ad-hoc checks. `src/rank_readers.py` and `src/figure2.py` read this file when `results/rt/` is absent.
- `screen/lists/`: the accession lists that define each screen set (RF1-less and RF2-less GTDB R232 bacteria at >= 70%
  CheckM2 completeness, genomes with stop-anticodon tRNAs, Pyl-machinery genomes, all Aristaeellaceae, Victivallales and
  Lentisphaeria, the equid and Tibetan-herbivore MAG projects, ...). Rerun any set with
  `uv run python src/batch_pooled.py <list>`.
- `screen/gtdb_r232_stops_release_factors.tsv.gz`: stop-codon usage and RF1/RF2/Pyl/Sec marker hits for all 199,923
  GTDB R232 species representatives (`src/gtdb_stream_scan.py`), the input to `src/analyze_gtdb_scan.py`.
- `screen/rank_readers_output.txt`: the ranking that found the readers. Saved before the two CNGB bins were added.
- `screen/site_concentration.tsv.gz`, `site_concentration_summary.txt`: the gene-specific recoding scan (Sec and Pyl
  recovered, no other candidate).

## Reference sequences
- `reference/rf.faa`: Swiss-Prot RF1/RF2 (for `rf_profile.py`). `reference/markers.faa.gz`: RF1, RF2, aRF1,
  Pyl/Sec machinery plus decoy families (for `gtdb_stream_scan.py`). `reference/arfb.faa`: ArfB. All are subsets of
  UniProtKB/Swiss-Prot.
