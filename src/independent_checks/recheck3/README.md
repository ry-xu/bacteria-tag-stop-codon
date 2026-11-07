Independent re-check of the three new Aristaeellaceae TAG readers by a separate agent (did not read src/, results/ or
RESULT.md), 2026-11-05 08:20-09:00. Tests A1 GCA_937907055.1, A2 GCA_045628485.1 (genus ZC3RG03), B GCA_030526225.1
(genus JAUNDF01); controls GCA_017577285.1, GCA_937886445.1, GCA_045617125.1 (ZC3RG03), GCA_934647165.1,
GCA_030529995.1 (CAKUFA01), E. coli K-12, B. subtilis 168.
Key outputs: out/summary_stops.tsv, out/rf_hmm.tsv, out/rf/prfA_synteny.tsv, out/tag_meaning_{A1,A2,B,ctrlA,ctrlB}.txt,
out/tag_swap_*.txt, out/trna/trna_report.txt. BLAST databases and out/cols (657 MB) not kept; regenerate with
scripts/run_all.sh <dir> (genomes from ../data/genomes/).
Result: TAG inside conserved alignments 7/17/15 vs TAA/TGA 0-3 (controls: 0 TAG); no RF1 in any test; prfA deleted
next to prmC with different junctions (A: prmC -> 170 nt -> glycoside hydrolase; B: prmC -> 30 nt -> sua5, control
has prmC-prfA-sua5); TAG -> Glu: Glu-Gln log-likelihood +64/+71/+102, mixture E 0.72-0.75; A1/A2 share a 76-nt CTA
tRNA (isotype Glu 65.1), none in B (bin incomplete: also lacks RF2, Tyr-GTA tRNA, 4 aaRS).
2026-11-05 10:20 (Ash): recheck4 found a bug in scripts/rf_fragments.py (locus filter used max() instead of min() on
e-values, dropping any locus that also had a weak hit). Fixed copy now in scripts/; rerun -> out/rf_fragments_fixed.txt:
A1/A2 carry only RF2 loci (plus a 38-aa, 37-bit RF1-query hit inside a YebA-like gene), B none, CAKUFA01 control RF2
only. Conclusions unchanged.
Note from the checker: no Pfam/TIGRFAM files (TIGR00019/TIGR00020, PF03462/PF00472) exist on this machine, so RF1/RF2
HMMs were built from 930 Swiss-Prot release factors (RF1 vs RF2 separated by >= 132 bits on Swiss-Prot), not checked
against the official models. The RF1 conclusions also rest on the deletion junctions and reciprocal Swiss-Prot hits.
