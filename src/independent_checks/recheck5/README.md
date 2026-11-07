Independent re-check of candidate lineage 6 (GCA_056209055.1, Tibetan wild ass, Aristaeellaceae; claimed UAG = Trp) by
a separate agent (did not read src/, results/ or RESULT.md), 2026-11-05 11:25-12:03; scoring code audited by a fresh
subagent (one bug fixed before the final numbers). Key outputs: out/tag_meaning5_T_coding.txt, out/tag_sites_T.tsv,
out/summary_stops.tsv, out/rf/rf_compare.txt, out/rf/prfA_synteny.tsv, out/rf/prfB_synteny.tsv,
out/trna/trna_trp_cta.txt, out/rt/*.readthrough.tsv. Regenerate with scripts/run_all.sh <dir> (~35 min, 3 threads).
Result: read-through TAG 27/47 vs TAA 0/719, TGA 1/890 (RF1-keeping relatives 0/19, 1/33; Glu reader 99/120). Only 3
TAGs inside conserved alignments (TAG is rare in this genome). TAG -> Trp on 21 sites where >= 2 controls agree (W 12,
Y 4, F 3, K 1, V 1; Q 0, E 0): W - Gln +98 [74,123], W - Glu +124, W - background +47, W - Tyr +37, W - "TAG kept at
former Trp sites without meaning" +17.6 [~2,32] (CI includes 0 with a looser site filter); Trp-vs-Tyr mixture 0.73.
TAGs sit at Trp columns as conserved as the genome's own TGG columns (e.g. GlmS HTRWATHG). RF1-like gene intact, GGQ,
accelerated (50-51% to relatives' RF1, 10-12th percentile of ~270 genes), motif VPAH-AGQK (PxT lost, 2-aa deletion).
RF2 not in the bin (assembly gap after rph). The "CUA" tRNA is tRNA-Trp-derived (86%) but tRNAscan-SE calls it
Undet (weak anticodon stem, 3 of 5 pairs): function unproven. Bin not clean (3.55 Mb vs 2.4-2.5; 4.5% duplicated;
42 conserved families missing incl. a ribosomal-protein operon), but no TAG evidence in the odd parts; not
contamination from the same-sample SFMI01 reader. GCA_934647165.1 lacks the whole prfA region (incomplete), so it is
not an RF1-keeping control. Verdict: moderately strong for Trp on few sites; decisive against Gln/Glu.
