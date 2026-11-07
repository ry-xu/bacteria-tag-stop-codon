Independent re-check of the Lachnospiraceae TAG readers (lineage 5) by a separate agent (did not read src/, results/ or
RESULT.md), 2026-11-05 09:53-10:14. Tests L1 GCA_056232285.1, L2 GCA_056261325.1 (Tibetan horse MAGs); controls
UBA1267 GCA_045640465.1, GCA_902803485.1, GCA_017959345.1; E. coli K-12, B. subtilis 168.
Key outputs: out/summary_stops.tsv, out/rf_hmm.tsv, out/rf/prfA_synteny.tsv, out/rf/rf_fragments_tests.txt,
out/tag_meaning_*.txt, out/tag_meaning4_*.txt, out/trna/trna_report.txt, out/trna/trna_compare.txt. Alignment tables
(out/cols) and BLAST/DIAMOND databases not kept; regenerate with scripts/run_all.sh <dir> (genomes: ../data/genomes/).
Result: TAG in conserved alignments 112 / 46 vs TAA/TGA 0-5 (controls 0-1 TAG); no RF1 (genes, HMM, raw-DNA tblastn);
controls keep bL31-yqhQ-prmC-prfA-mgsA; in L2 prmC -> CUA tRNA -> mgsA (the tRNA sits where prfA was); L1's join is
split over contigs. TAG -> Gln: unanimous-control positions Q 72% / 69% (CAG -> Q 81% / 79%, GAG -> E 88%), Gln-Glu
log-likelihood +469 / +190, mixture Gln 0.90 / 0.89. CUA tRNAs 83% identical between L1 and L2, isotype Gln/Ile2/Cys
(not Glu); GlnRS present, no GatCAB. L1 and L2 share 16 TAG positions (3.5x chance): one origin. Caveats: low-level
contamination in L2 (a divergent rplK-rplA copy, no TAGs); L2 uses TAG 2.3x less than L1. Found and fixed a bug in
recheck3's rf_fragments.py (see ../recheck3/README.md).
