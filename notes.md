# Notes: bacteria that read the TAG stop codon as an amino acid

The agents' full write-up (method, evidence, checks, open points). The short version is in [`README.md`](README.md).

**Summary.** We screened about 41,000 bacterial genomes for genes whose protein continues past the stop codon.
Twelve genomes read through TAG in 38 to 347 genes each, and through TAA or TGA in only 0 to 7. They fall in at least
six independent lineages: one in Victivallales (Verrucomicrobiota), four genera of Aristaeellaceae and one
Lachnospiraceae (both Bacillota, class Clostridia). All twelve lack a working release factor 1 (RF1, the protein that
ends translation at UAG and UAA). Nine have lost the prfA gene; wherever the locus is assembled, each lineage has a
different deletion breakpoint. The other three keep an RF1 whose stop-recognition motif is broken. Their closest relatives keep RF1 and use TAG as an
ordinary stop. From the residues that relatives carry at TAG positions, UAG appears to encode glutamate in three
Aristaeellaceae genera and glutamine in the Lachnospiraceae lineage. In a fourth Aristaeellaceae genus (HGM11604, the
weakest case) it probably encodes tryptophan. In Victivallales its meaning is unresolved. Ten of the twelve genomes
come from the guts of hindgut-fermenting mammals (horses, Tibetan wild asses, a warthog). **All evidence is
computational, from uncultured metagenome-assembled genomes (MAGs). No protein data exist for any of these organisms.**
We are posting it so that experts can check it.

## What was known
- **Bacteria.** Every bacterial stop-codon reassignment reported so far is at UGA: UGA to Trp (NCBI code 4:
  Mycoplasmatales and Spiroplasma; also some Patescibacteriota and, in 2026, Actinomycetota, Parks et al. 2026) and
  UGA to Gly (code 25: Gracilibacteria/Absconditabacteria). UAG is read as pyrrolysine in a few bacterial genes, but
  only in those genes. RF1-deleted *E. coli* exists only in the lab.
- **TAG read genome-wide.** Known in eukaryotes: Gln in ciliates, diplomonads and dasycladacean green algae; Glu in
  some ciliates and *Blastocrithidia*; Tyr in *Mesodinium*. Also known in the plastids of the parasitic plant
  *Balanophora* (Trp, code 32), in gut phages (mostly Gln), and in archaea. Kivenson et al. (*Science* 2025) described
  an archaeal code in which every TAG encodes pyrrolysine (code 34).
- **No report found.** We found no report of a natural bacterium that reads TAG as an amino acid genome-wide, or that
  has lost RF1. Searches up to 2026-11-06 covered the web, Europe PMC (including bioRxiv), the 2026 code-prediction
  tools (gTranslate, KACI, Codetta) and gTranslate's published per-genome code calls, which list four of these
  genomes as standard code 11. A paper in press or a talk could still have been missed.

## The six lineages and twelve genomes
Read-through = genes read through / genes tested, from the full (unsampled) test. Completeness/contamination is
CheckM2 from GTDB R232 unless marked. Full table with more columns: [`data/genomes.tsv`](data/genomes.tsv).

| Lineage | Genome | Host, place | BioProject | Compl./contam. % | TAG / TAA / TGA read-through | RF1 | UAG meaning (confidence) |
|---|---|---|---|---|---|---|---|
| L1 Victivallales (Verrucomicrobiota), fam. UBA1829, gen. CASERY01 | GCA_949290385.1 | chicken caeca (49-sample co-assembly), Ethiopia | PRJEB57055 | 78.5 / 2.0 | 281/739 · 4/496 · 5/509 | lost (bL31 to ywjA, 402 bp) | unresolved |
| L1 | GCA_965332555.1 | red kangaroo faeces, Australia, 2023 | PRJEB89408 | 86 / 1 (depositor) | 181/484 · 3/435 · 0/376 | lost (locus rearranged) | unresolved |
| L2 Aristaeellaceae SFMI01 | GCA_024404275.1 | horse faeces (11-sample co-assembly), France, 2015 | PRJNA438436 | 74.7 / 3.6 | 92/181 · 3/445 · 2/331 | lost (lepA to topIII, ~20 bp) | Glu (strong) |
| L2 | GCA_045599485.1 | Tibetan wild ass faeces, Tibet, 2019 | PRJNA1030165 | 81.7 / 1.7 | 201/376 · 3/566 · 3/492 | lost (same junction) | Glu (strong) |
| L3 Aristaeellaceae ZC3RG03 | GCA_937907055.1 | horse gut, Inner Mongolia, 2013 | PRJEB50943 | 89.0 / 4.7 | 45/130 · 5/510 · 5/676 | lost (prmC to glycoside hydrolase, 170 nt) | Glu (moderate) |
| L3 | GCA_045628485.1 | Tibetan horse faeces, Tibet, 2019 | PRJNA1030165 | 89.8 / 4.4 | 58/131 · 6/489 · 6/676 | lost (same junction, independent assembly) | Glu (moderate) |
| L4 Aristaeellaceae JAUNDF01 | GCA_030526225.1 | warthog faeces, Zoo Saarbrücken, 2022 | PRJNA983076 | 79.8 / 2.0 | 106/260 · 2/698 · 4/356 | lost (prmC to sua5, 30 nt) | Glu (strong; decoding tRNA not in bin) |
| L5 Lachnospiraceae near UBA1267 | GCA_056232285.1 | Tibetan horse faeces, Tibet, 2019 | PRJNA1030165 | not estimated | 347/690 · 4/512 · 1/184 | lost | Gln (strong) |
| L5 | GCA_056261325.1 | Tibetan horse faeces, Tibet, 2019 | PRJNA1030165 | not estimated | 167/405 · 3/532 · 0/280 | lost (prmC to CUA tRNA to mgsA) | Gln (strong) |
| L6 Aristaeellaceae HGM11604 | GCA_056209055.1 | Tibetan wild ass faeces, Tibet, 2019 | PRJNA1030165 | not estimated (oversized bin) | 42/154 · 3/795 · 7/891 | present, PxT motif broken (PAH) | probably Trp (weak) |
| L6 | QTP_TA20190011_bin.237 | second wild ass, same site | CNGB CNP0001390 | not estimated | 38/133 · 3/671 · 6/785 | same broken RF1 | probably Trp (weak) |
| L6 | QTP_TA20190013_bin.296 | third wild ass, same site | CNGB CNP0001390 | not estimated (chimeric bin) | 43/158 · 7/890 · 7/949 | same broken RF1 | probably Trp (weak) |

L1 is one lineage seen twice (68.8% AAI apart). L2 and L5 are each two species of one lineage. L3 is one species from
two horse populations six years apart. L6 is one species seen in three animals: the three genomes are 99.5-99.6%
identical, so they are not three independent observations. The two L6 bins are not in NCBI. They come from the
depositors' full MAG set on CNGB (Li et al., *Microbiome* 2025, 13:236).

## Caveats (read these first)
- **Computational only.** No isolate, no culture, no proteomics, no ribosome profiling. Every amino-acid assignment is
  an inference from sequence.
- **Uncultured MAGs.** All twelve genomes are metagenome bins: incomplete (75-90% where estimated), with some
  contamination. RF1 absence from a bin means little on its own. The evidence for loss is the deletion junction inside
  a contig, which we have for L1, L2, L3, L4 and one L5 genome. JAUNDF01 (L4) is missing RF2 from its bin as well. For
  L5 and L6 no completeness estimate exists, and two L6 bins are oversized or chimeric.
- **Meanings are inferences, and they differ in strength.** Glu (L2, L4) and Gln (L5) fit nearly as well as a real
  codon for those amino acids does in the same test. L3's Glu call rests on fewer sites: 8 of 17 at conserved
  positions. L6's Trp call is weak (see below). L1 does not fit any single amino acid.
- **Lineage 6 is the weakest of the six on every axis.**
  - Its read-through signal is the lowest.
  - RF1 is still present. Only its PxT stop-recognition motif is broken, and nothing shows the factor has stopped
    acting on UAG while still acting on UAA, so UAG may be read ambiguously.
  - Its CUA tRNA is called "Undet" by tRNAscan-SE, so its function is unproven.
  - Only ~55-60% of TAG sites look Trp-like, against 82-95% Glu-like for the best Glu reader.
  - The real competitor to Trp is "TAG sits at former Trp positions without meaning Trp". Trp beats it by only
    +13 log-likelihood units, 95% interval [-3, 24], which includes zero.
  - Seven fixed TAGs sit at Tyr positions.
- **Read counts in the L6 raw-read check are inflated.** Mates were mapped as unpaired reads, without duplicate
  removal. The conclusion that the TAGs are fixed in four animals stands.
- **Screen gaps.**
  - About 100 genomes across the screen failed to download (7 of them in the RF1-less set), so "all 9,522 RF1-less
    GTDB genomes" means all but a handful.
  - Most genomes were tested on a deterministic subsample of genes.
  - One L6 relative, GCA_056272465.1 (Tibetan antelope, RF1 not in bin, 5 of 19 TAG genes read through on the
    subsample), was never given a full run.
  - No shared-TAG-position test was run between L6 and the other lineages. They are 54-63% AAI apart and differ in
    mechanism and tRNA.
- **Phylogeny is GTDB's.** We did not build our own trees.
- **Literature can be missed.** See "What was known".

## The evidence
**1. Read-through is TAG-specific** ([`src/readthrough.py`](src/readthrough.py)).
- Method: genes are called with pyrodigal (code 11). Each gene is extended past its stop to the next in-frame stop
  (up to 300 codons), and the joined protein is aligned to Swiss-Prot with DIAMOND. A gene reads through if one
  alignment covers its stop with at least 15 residues on each side.
- The twelve read through TAG in 38-347 genes (27-53% of TAG-ending genes) and through TAA or TGA in 0-7 genes.
- RF1-keeping relatives read through TAG in 0-9 genes ([`data/readthrough_counts.tsv`](data/readthrough_counts.tsv)).
- The read-through genes are core genes: EF-G, RNA polymerase, gyrase, aminoacyl-tRNA synthetases, ribosomal
  proteins, and RF2 itself. In-frame TAGs are spread over many contigs.
- Ranking ([`src/rank_readers.py`](src/rank_readers.py)) requires both a high rate and a high absolute,
  sampling-corrected excess over the other two stops. Ranking by rate alone favours noisy MAGs that rarely use TAG.

**2. RF1 is gone, at a different breakpoint in each lineage** ([`data/rf1/breakpoints.tsv`](data/rf1/breakpoints.tsv)).
Each RF1-keeping relative has prfA at the position marked in the junction.
- L1: bL31 is followed directly by ywjA (402 bp).
- L2: lepA is followed by topIII (~20 bp), in both genomes.
- L3: prmC runs 170 nt into a glycoside hydrolase. The same junction appears in two independent assemblies, one on an
  83.6 kb contig.
- L4: prmC is followed by sua5 after 30 nt, inside a 38 kb contig.
- L5: in GCA_056261325.1, prmC is followed by the CUA tRNA and then mgsA, so the UAG-reading tRNA sits where prfA was.
- L6: RF1 is still at its native locus, intact with GGQ, but only 43% identical to Swiss-Prot RF1 (relatives 54-58%).
  Its PxT motif reads PAH.
- RF2 (SPF motif) is present in every genome whose RF2 region was assembled. One L3 genome reads PPF for SPF at one
  position, likely an assembly or population variant.
- Lineages L2, L3, L4 and L6 are four genera of one family, each nested among RF1 keepers. L2, L3 and L4 are
  0.34-0.45 substitutions/site apart in the GTDB R232 tree and share 0 of 41 TAG positions, against 40 of 48 between
  the two L2 genomes. A single ancestral loss would need RF1 to be regained, at the same locus, in hundreds of
  genomes.

**3. The TAGs are in the raw reads** ([`data/read_support/`](data/read_support)). Sites are in-frame TAGs covered by
at least 3 reads; "fixed" means TAG in at least 90% of reads.

| Lineage | Sample | Sites covered | Fixed | Mixed | Absent |
|---|---|---|---|---|---|
| L1 | kangaroo, ERR14999989 | 583 of 730 | 563 | 20 | 0 |
| L2 | one of the 11 pooled horses, SRR17543912 | 122 of 248 | 91 | 19 | 12 mostly another codon* |
| L3 | horse, SRR6425796 | 199 of 201 | 179 | 20 | 0 |
| L4 | warthog, SRR24915835 | 355 of 359 | 340 | 14 | 1 |
| L6 | four wild asses, CNCB GSA CRA013436 | 209 of 213 | 158 fixed in every animal that covers them | | none, anywhere |

\* This horse is one of 11 pooled in the co-assembly. At non-fixed L2 sites the other codon is a Glu codon in 66 of
93 reads (GAG 54, GAA 12). This fits Glu but adds little to the ortholog test: at a polymorphic site the other
codon is simply what relatives carry there, whatever TAG means (recheck6's critique, which led us to withdraw the
same argument for L6).

L5's raw reads were not found in SRA/ENA. The TAGs are not assembly errors.

**4. Independent re-checks.** Six separate checks were each run by an agent that had not seen the main pipeline's code
or results ([`data/independent_checks/`](data/independent_checks), [`src/independent_checks/`](src/independent_checks)).
They used tblastn of ~456 conserved *E. coli* / *B. subtilis* proteins (or of relatives' proteomes, for L6) against raw
genome DNA, with no gene calls. Stops counted inside alignments:

| Genome | TAG | TAA / TGA |
|---|---|---|
| L1 chicken / kangaroo | 51 / 37 | 4/1, 3/4 |
| L2 horse / wild ass | 19 / 48 | 2/0, 0/1 |
| L3 horse / Tibetan horse | 7 / 17 | 3/2, 3/3 |
| L4 warthog | 15 | 2/0 |
| L5 two Tibetan horses | 112 / 46 | 5/2, 4/0 |
| L6 three genomes (relatives' proteomes) | 51 / 43 / 51 | 0-3 / 2-4 |
| 14 RF1-keeping closest relatives (L1-L5) | 0-1 each | |

The method recovers *E. coli*'s selenocysteine TGAs. *B. subtilis* gives 3 TAGs, at alignment edges. Further checks:
- A homology-free test for L6: 76-78% of TAG-ending genes keep an open frame for at least 50 more codons, against
  23-36% for TAA/TGA-ending genes.
- TAGs sit on many contigs, with the bin's GC and coverage.
- No check found signs of contamination, pseudogenes, selenocysteine, pyrrolysine (no PylS/PylBCD) or a
  minority-suppressor explanation.

**5. What TAG encodes.**
- Method ([`src/ortholog_profiles.py`](src/ortholog_profiles.py)): for every codon of the reader, take the residue at
  the orthologous position in close RF1-keeping relatives. Then ask which amino acid's own codons give a residue
  profile most like TAG's, by log-likelihood and a mixture fit. Real codons give the calibration: a real codon scores
  40-67% "its own" amino acid in the same test.
- Per-site tables: [`data/tag_sites/`](data/tag_sites). Outputs: [`data/meaning/`](data/meaning).

| Lineage | TAG sites | Relatives' residue at TAG | Calibration (real codon) | Best vs next (log-likelihood) | Independent re-check |
|---|---|---|---|---|---|
| L2 | 169 / 72 | Glu 49.7% / 52.8% (background 6.3-6.4%) | ~50% | Glu beats Asp by 149 / 68; mixture Glu 0.85 | Glu 71% in both (GAG gives 72-73%). Conserved-protein majority: 34 of 47 sites Glu |
| L3 | 44 / 46 | Glu 45.5% / 41.3% (background 7.0%) | 44-63% | Glu beats Gln by 28 / 23; mixture 0.68 / 0.61 | Glu 52-57%, Gln 6-8%. Only 8 of 17 conserved sites Glu |
| L4 | 89 | Glu 47.2% (background 6.6%) | 44-63% | Glu beats Gln by 73; mixture 0.76 | 10 of 13 conserved sites Glu |
| L5 | 367 / 140 | Gln 31.1% / 40.0% (background 2.8%) | CAG gives Gln 41% / 40% | Gln beats Glu by 242 / 122; mixture 0.73 / 0.89 | Gln 71.7% / 69.0% where all controls agree (CAG gives 81% / 79%) |
| L6 | 52 | Trp 48% | TGG gives Trp 67% | Trp beats Tyr by 96, Gln by 115; mixture Trp 0.65, Tyr 0.08 | Trp at 26/51, 23/43, 27/51 sites (TGG gives 85-86%): ~55-60% Trp-like |
| L1 | 289 | Ala 19.7% (2.1x), Gln 8.3% (3.1x) | | background beats every single amino acid; mixture Ala 0.29, Gln 0.18 | unresolved |

- **Could TAG just have arisen from a Glu codon by mutation?** That would place TAG at Glu positions without meaning
  Glu. Three things argue against it:
  - At L4's TAG sites the relatives' codon is GAA 13 times and GAG 8 times. GAA is two changes from TAG.
  - CAG and AAG are equally one step from TAG, yet Gln and Lys sit at background.
  - Between the two L2 genomes, TAG swaps with GAG/GAA at 14 of 37 positions and with CAG at 1.
- **The reassignments are old.**
  - The two L1 genomes share 33 TAG positions against 1.3 expected by chance.
  - The two L5 species share 16, about 3.5 times chance.

**6. tRNAs** ([`data/trnas/`](data/trnas)).
- L2: a 76-nt tRNA that pairs with UAG and most resembles a glutamate tRNA (tRNAscan-SE isotype scores Glu 70.3,
  Gln 59.0). It is identical in both genomes and absent from relatives.
- L3: a different 76-nt CUA tRNA at a different locus, also Glu-like (Glu 65.1, Cys 58.6). Unrelated to L2's.
- L5: a CUA tRNA most like Gln/Ile2/Cys tRNAs, sitting where prfA was. These genomes have GlnRS and no GatCAB, which
  fits glutamine.
- L6: a CUA copy of the genome's own tRNA-Trp, 2 kb away. tRNAscan-SE calls it "Undet".
- L1: no CUA tRNA in either bin. 125 M kangaroo-sample reads contain CUA tRNAs, but none is linked to the MAG or
  resembles the clade's own tRNAs.
- L4: no CUA tRNA in the 80%-complete bin.
- Charging of these tRNAs is unproven. Bacterial GluRS reads anticodon base 36, which is A in a CUA anticodon.

**7. The screen** ([`data/screen/`](data/screen)).
- 42,230 genome accessions in all: about 41,200 bacterial genomes plus ~900 archaeal controls. They include:
  - all 9,522 RF1-less GTDB R232 bacteria at 70% or more completeness;
  - 9,489 RF2-less bacteria;
  - 6,713 genomes with stop-anticodon tRNAs;
  - 919 genomes with pyrrolysine machinery;
  - all 1,705 Aristaeellaceae and 1,228 Lentisphaeria;
  - 2,149 equid-gut MAGs;
  - 12,411 other Tibetan-herbivore MAGs.
- The screen recovered known codes blind: code-34 archaea; UGA = Trp in *Equivita*/*Gorillivita* (Parks et al. 2026)
  and in Minisyncoccia; UGA = Gly in Gracilibacteria.
- No other bacterial genome reads TAG genome-wide.
- A side scan for gene-specific recoding that recurs across phyla found only selenocysteine (TGA) and pyrrolysine
  (TAG): no candidate 23rd amino acid.

## What would confirm it
- **Proteomics on the original gut samples.** Peptides spanning TAG sites in highly expressed proteins (EF-Tu, EF-G,
  ribosomal proteins) would settle both the read-through and the amino acid. The best targets are the Tibetan wild-ass
  and horse samples of PRJNA1030165 (L2, L3, L5, L6) and the warthog sample (L4). No suitable public metaproteome
  exists. The only chicken-caecum metaproteome on PRIDE (PXD044520) is from gnotobiotic chicks.
- **Culturing.** Any of these lineages would allow tRNA charging assays, ribosome profiling and reporter tests. L2's
  SFMI01 occurs in both horses and wild asses.
- **Long-read assembly** of the same samples, to close the prfA loci and check bin completeness, especially for L5 and
  L6.
- **For L1**, finding the UAG-decoding tRNA.

## How to reproduce
Requirements:
- macOS or Linux, [uv](https://docs.astral.sh/uv/), and `diamond` on PATH (`brew install diamond`).
- For the raw-read checks: `minimap2` and `sratoolkit`.
- For the independent checks: BLAST+, and tRNAscan-SE 2 (set `TRNASCAN_BIN`).

```bash
uv sync                                   # Python deps (pyrodigal, pyaragorn, requests, matplotlib, biopython, pillow)
scripts/setup_reference_data.sh           # Swiss-Prot -> DIAMOND db, GTDB R232 metadata + trees (~0.6 GB)
uv run python src/readthrough.py GCA_030526225.1   # one genome, ~1 min: expect TAG 106/260, TAA 2/698, TGA 4/356
scripts/fetch_cngb_bins.sh                # the two L6 bins from CNGB (slow server)
scripts/reproduce_core.sh                 # read-through, RF1, ortholog tests for all 12; screen ranking; figure
```

- `GENOME_DIR` and `SPROT_DB` point the code at an existing genome cache or DIAMOND database.
- Every script in `src/` has a usage line at the top. Run them as `uv run python src/<script>.py ...`.
- The full screen is `src/batch_pooled.py` over the lists in `data/screen/lists/`. It takes days on one machine and
  downloads each genome once. `src/gtdb_stream_scan.py` streams 187 GB of GTDB gene files without storing them.

**Checked while building this repository** (2026-11-06, against the original outputs):
- `readthrough.py` on JAUNDF01 gave identical counts.
- `ortholog_profiles.py` and `ortholog_codons.py` on JAUNDF01 gave byte-identical output.
- `ortholog_profiles.py` on SFMI01 and Victivallales reproduced the numbers above.
- `rank_readers.py` on the shipped screen table gave the same reader list. The only differences were the two later
  CNGB bins and the order of tied genomes.
- `figure2.py` redrew the figure (40,998 other genomes, 100 known UGA readers).

**Not reproducible from this repository alone:**
- `src/extract_stop_trnas.py` needs a non-public archived ARAGORN screen of all GTDB R232 genomes. The genome list it
  produced is shipped (`data/screen/lists/armB_genomes.tsv`).
- `site_clusters.py`, `site_concentration.py` and `site_detail.py` need the per-gene screen outputs (`results/rt/*.tsv`,
  337 MB in 42k files). These are not shipped; rerun the screen to regenerate them. Their summaries are shipped.
- `sample_asm_scan.py` needs the 49 per-sample assemblies behind CASERY01 (ENA, PRJEB57055).
- The independent checks (`src/independent_checks/`) have their paths parameterized but were not rerun here.
- Swiss-Prot changes between releases. The original run used the release current on 2026-11-04, so individual hits
  can differ slightly.

## Repository layout
- `src/`: analysis code. Ash's pipeline, lightly adapted: paths made relative, genome cache and Swiss-Prot path
  overridable, `rank_readers`/`figure2` able to read the shipped screen table. Also contains `independent_checks/` (the
  re-check scripts and their shared inputs) and `trna_hunt/`.
- `data/`: small result tables; see [`data/README.md`](data/README.md).
- `scripts/`: setup, CNGB download, core reproduction.

## Credit
This work was done by AI research agents in Ryan Xu's lab, and the analysis was done by Claude agents (Anthropic's
Claude models running as autonomous research agents):
- An agent called Ash ran the screen and the main analysis, and found the lineages.
- Separate Claude agents ran the six independent re-checks and two pre-publication fact-checks (prior art, numbers,
  wording).
- This repository was assembled by a Claude agent.

Ryan Xu is the human responsible for this work. No outside expert has reviewed it yet.

Data come from public deposits. We thank the depositors of these MAGs and read sets:
- PRJEB57055 (Addis Ababa University);
- PRJEB89408 (University of Queensland);
- PRJNA438436;
- PRJNA1030165 and CNP0001390 (Li et al. 2025);
- PRJNA983076;
- PRJEB50943;
- and GTDB, UniProt and NCBI.

**Contact:** ryan@ryanxu.net.

## License
- Code: MIT; text, data tables and figures: CC BY 4.0 (see [`LICENSE`](LICENSE)).
- Third-party material keeps its own terms. Files in `data/reference/` are subsets of UniProtKB/Swiss-Prot (CC BY 4.0).
  The GTDB-derived tables in `data/screen/` are computed from GTDB R232; check GTDB's terms before redistributing them.

## References
- Kivenson V, Peters SL, Borrel G, Kivenson A, Roe LT, Hamlish NX, Fadhlaoui K, Schepartz A, Gribaldo S, Hettich RL,
  Banfield JF. An archaeal genetic code with all TAG codons as pyrrolysine. *Science* 390 (2025).
  doi:10.1126/science.adu2404
- Parks DH, Chaumeil P-A, Chuvochina M, Hugenholtz P. Stop codon reassignment to tryptophan in members of the bacterial
  phylum Actinomycetota. *Microbial Genomics* 12 (2026). doi:10.1099/mgen.0.001767
- Chaumeil P-A, Hugenholtz P, Parks DH. gTranslate. *PeerJ* (2026). doi:10.7717/peerj.21825
- Fakih et al. tRNAs with 4-bp anticodon stems. *Nucleic Acids Research* (2026). doi:10.1093/nar/gkag327
- Ivanova NN et al. Stop codon reassignments in the wild. *Science* 344 (2014).
- Su H-J et al. Novel genetic code and record-setting AT-richness in the highly reduced plastid genome of the
  holoparasitic plant *Balanophora*. *PNAS* 116 (2019).
- Mukai T et al. Highly reproductive *Escherichia coli* cells with no specific assignment to the UAG codon.
  *Scientific Reports* 5 (2015).
- Kachale A et al. Short tRNA anticodon stem and mutant eRF1 allow stop codon reassignment. *Nature* 613 (2023).
- Li et al. A unified catalog of microbial species reference genomes from the Qinghai-Tibet Plateau. *Microbiome* 13,
  236 (2025).
- Shulgina Y, Eddy SR. A computational screen for alternative genetic codes in over 250,000 genomes. *eLife* 10 (2021).
