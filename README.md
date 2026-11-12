> **This repository is frozen.** The work now lives on, together with later results, at
> [github.com/ry-xu/genome-research](https://github.com/ry-xu/genome-research) (folder `tag-bacteria`).

# Bacteria that read the stop codon TAG as an amino acid

Almost every gene ends in one of three stop codons: TAA, TAG or TGA. Bacteria were known to rewrite TGA into an amino acid
(tryptophan in *Mycoplasma*, glycine in Gracilibacteria), but never TAG.

We screened about 41,000 bacterial genomes and found 12, in six independent lineages, that read TAG as an amino acid across
the whole genome. Five of the six come from the guts of horses, Tibetan wild asses and a warthog. In three of them, TAG most
likely means glutamate.

---

- **AI-generated.** The screen, the analysis and the write-up were done by Claude (Anthropic) agents, run and curated by Ryan Xu.
- **Checked by AI, not yet by a human expert.** A second agent re-checked every lineage with a different method (tblastn
  against the raw genome DNA, no gene calls), and two more agents fact-checked the claims and searched for prior work.
- **Computational only.** The genomes were assembled from gut DNA by other groups (metagenome-assembled genomes, MAGs). None
  of these bacteria has been grown, and there is no protein data for them. The amino acid assignments are inferences.
- Corrections and criticism are very welcome: ryan@ryanxu.net.

### What was found?

| lineage | genomes | from | genes through TAG | through TAA / TGA | RF1 | TAG means |
|---|---|---|---|---|---|---|
| Aristaeellaceae, genus SFMI01 | 2 | horse (France), Tibetan wild ass | 92/181, 201/376 | 2–3 | deleted | glutamate |
| Aristaeellaceae, genus ZC3RG03 | 2 | horses (Inner Mongolia, Tibet) | 45/130, 58/131 | 5–6 | deleted | glutamate |
| Aristaeellaceae, genus JAUNDF01 | 1 | warthog (zoo, Germany) | 106/260 | 2–4 | deleted | glutamate |
| Aristaeellaceae, genus HGM11604 | 3 | Tibetan wild asses | 38–43 of 133–158 | 3–7 | present but broken | probably tryptophan |
| Lachnospiraceae, near UBA1267 | 2 | Tibetan horses | 347/690, 167/405 | 0–4 | deleted | glutamine |
| Victivallales, family UBA1829 | 2 | chicken (Ethiopia), red kangaroo (Australia) | 281/739, 181/484 | 0–5 | deleted | not yet known |

"Genes through TAG" counts, of the genes that end in TAG, how many make a protein that carries on past it and keeps matching
its relatives. `notes.md` numbers the lineages L1–L6 in the order they were found (SFMI01 = L2, ZC3RG03 = L3, JAUNDF01 = L4,
HGM11604 = L6, Lachnospiraceae = L5, Victivallales = L1).

### How do we know TAG isn't a stop here?

- **It's only TAG.** In these genomes, 38 to 347 genes run straight through TAG, essential genes included (DnaA, RNA
  polymerase, EF-G). Through TAA or TGA it's 0 to 7, and in their closest relatives about zero.
- **The stop protein is gone.** All twelve lack a working release factor 1 (RF1, the protein that ends translation at TAG).
  Nine have deleted the gene; the three HGM11604 genomes keep a copy with its stop-recognition motif broken. Their closest
  relatives keep RF1 and use TAG as an ordinary stop.
- **Six separate times.** Each lineage sits among relatives that still stop at TAG, and each lost or broke RF1 in its own way
  (different deletion breakpoints; in HGM11604 a copy that decayed in place).
- **The TAGs are real.** Raw sequencing reads support them, and public horse-gut RNA data show the genes transcribed with
  their TAGs.

### Which amino acid does TAG code for?

For each TAG inside a gene, we look at the same position in the same protein in close relatives that still stop at TAG,
and note the amino acid they have there. Proteins drift, so even a real codon doesn't match its relatives every time: in this
test a real glutamate codon shows glutamate at about 45–65% of its sites, and a random site shows it at 6–7%.

- **Glutamate** (SFMI01, ZC3RG03, JAUNDF01): relatives have glutamate at 41–53% of TAG sites, against 6–7% at random sites,
  the same range as a real glutamate codon. It beats every other amino acid clearly, and the independent re-check agrees.
- **Glutamine** (Lachnospiraceae): glutamine at 31–40% of TAG sites, against 3% at random; at the sites where all controls
  agree, 69–72%, against 79–81% for a real glutamine codon.
- **Probably tryptophan** (HGM11604): tryptophan at about half the sites (26 of 51), where a real tryptophan codon gives 85%.
  A lean, not a result.
- **Unresolved** (Victivallales): no amino acid stands out.
- **Could TAG have simply come from a glutamate codon by mutation?** GAG → TAG is one letter, so TAG could sit at glutamate
  positions without meaning glutamate. Against this: at JAUNDF01's TAG sites the relatives' codon is more often GAA (13 of
  21), two letters from TAG; and CAG and AAG, equally one letter from TAG, don't show glutamine or lysine above background.
- **tRNAs.** SFMI01 and ZC3RG03 each carry a new tRNA that pairs with TAG and is built like a glutamate tRNA. The
  Lachnospiraceae lineage and HGM11604 each have their own TAG-pairing tRNA (HGM11604's is derived from a tryptophan tRNA).
  None has been found yet in JAUNDF01 (an incomplete genome) or Victivallales.

### What was known?

TAG is read as an amino acid across whole genomes in some single-celled eukaryotes (glutamine or glutamate), and in archaea
that read every TAG as pyrrolysine, a rare 22nd amino acid ([Kivenson, Banfield et al., *Science* 2025](https://doi.org/10.1126/science.adu2404)).
Earlier screens of bacteria ([Shulgina & Eddy 2021](https://elifesciences.org/articles/71402); Parks et al. 2026) found
stop codon changes only at TGA. Some bacteria use TAG for pyrrolysine, but only in a handful of genes.

### What would settle it?

Protein data. The archaeal code was confirmed by finding pyrrolysine in measured proteins; finding glutamate at TAG positions
in proteins from horse-gut samples would do the same here. No public metaproteome from equid guts exists yet.

### Files

- `notes.md`: the agents' full write-up (method, evidence, controls, checks, open points, references)
- `data/`: the 12 genomes (`genomes.tsv`), read-through counts, every TAG site with its relatives' residues, RF1 loci, tRNAs,
  read support, the independent re-checks, and the full screen; see [`data/README.md`](data/README.md)
- `src/`, `scripts/`: the analysis, runnable with `uv run`; `scripts/` downloads the large reference files
