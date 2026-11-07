# recheck6: skeptic re-test of candidate lineage 6 (UAG read as sense, probably Trp)

This is an independent re-test done on 2026-11-05. It covers GCA_056209055.1 and the two new CNGB bins, QTP_TA20190011_bin.237 and
QTP_TA20190013_bin.296. I used my own methods and did not edit any existing files.

## Verdict
- **Read-through genome-wide: HOLDS.**
  - **Same species.** The bins are the candidate's species: 99.5-99.6% identity over 93.5% (bin.237) and 85.3% (bin.296) of their length. That is closer than the claimed ~97%.
  - **TAGs inside alignments.** I ran tblastn of three relatives' proteomes against raw DNA and counted stops inside high-identity alignments (HSP >= 50% id, both 30-aa flanks >= 40% id). TAG/TAA/TGA counts:
    - candidate 51/1/3, bin.237 43/0/2, bin.296 51/3/4;
    - RF1-keepers 3/2/1 and 0/1/3;
    - GCA_056207115.1 5/1/0;
    - Glu reader 112/0/1.
  - **Homology-free test.** I measured how far the reading frame stays open after the stop of genes ending in TAG. In the three bins, 76-78% have >= 50 open codons after the TAG, against 23-36% after TAA/TGA. The controls are 30-48%, and the Glu reader is 83%.
  - **Sites shared across assemblies.** Of the candidate's TAG sites found in each new bin, 161/163 (bin.237) and 200/201 (bin.296) are also TAG there.
- **Same decayed RF1 and tRNA, and the contig_000065 chimera: HOLDS.**
  - The decayed RF1 is identical at the nucleotide level in all three bins. Where relatives have PxT it reads VPA--H; GGQ is kept.
  - All three carry the same 72-nt tRNA with a CTA anticodon, 2.1 kb from tRNA-Trp. tRNAscan-SE calls it "Undet" (43 bits), so its function is unproven.
  - Positions 72.5 kb to the end of contig_000065, which include the normal RF1 (PIT motif), are 99.3% identical to TA20190008_bin.73. That is a genome of GCA_056207115.1's species: 99.75% ANI.
  - bin.73's RF1 is 100% identical to the normal RF1 and fills the gap in GCA_056207115.1 (= TA20190013_bin.200, the same animal as bin.296).
  - bin.296's own prfA locus, on contig_000213, carries the decayed RF1.
  - Note: GCA_056207115.1 has two high-scoring Sup-CTA tRNAs (64 bits) yet does not read through.
- **Read support: WEAKENED (the interpretation, not the data).**
  - The coordinates and strands are correct: all 815 site-rows read TAG in the reference.
  - Flaw: the f1 and r2 reads are mapped as unpaired, independent reads, and duplicates are not removed, so read counts are inflated.
  - Flaw: the NM filter ignores soft-clipped bases. This is harmless here: no other species in the four animals is >= 95% identical around any polymorphic site (`out/crossmap_summary.txt`).
  - The TGG "control" mixed rate is just sequencing error.
  - Over-interpretation: at polymorphic sites the alternative codon is the one relatives carry at that column (TGG at 4/4 Trp columns; CAG at Gln, Met and Asp columns). That is the ancestral allele and is expected under any meaning. The absence of TAA says nothing about meaning.
  - By the claimant's own Glu-reader logic, CAG (10) >= TGG (9) would point to Gln as much as to Trp.
  - What stands: TAG is fixed at most sites in four animals.
- **Meaning: WEAKENED to "probably Trp".**
  - At TAG sites the relatives' residue is W in 26/51 (candidate), 23/43 (bin.237) and 27/51 (bin.296); Q in 1-3 and Y in 4-7.
  - Real TGG codons give W in 85-86% of the same alignments. So TAG sites look only about 55-60% Trp-like. The Glu reader's TAG sites look 82-95% Glu-like.
  - TAG columns are as W-conserved in Swiss-Prot as the genome's own TGG columns: mean W fraction 0.27 vs 0.28, from 16 informative sites.
  - Hand spot-checks find real conserved Trp: TrmB GKW, OmpR-family DRVW (W in 48/50 Swiss-Prot hits), β-galactosidase LDW. They also find 7 fixed TAGs at Tyr columns.

## Rerun
Run `scripts/run_all.sh`. It takes about 15 minutes, needs less than 2 GB of RAM, and rebuilds and then deletes the temporary BLAST database of 584 CNGB bins. For the chimera check, blastn `out/rf/c65_65k.fna` (positions 65 kb to the end of contig_000065) against the candidate, bin.237 and CNGB TA20190008_bin.73.

Key outputs:
- `out/stops/summary.txt`
- `out/downstream_orf.tsv`
- `out/*contigs.tsv`
- `out/rf/rf_hits.tsv`
- `out/trna/`
- `out/read_sites_check.txt`
- `out/crossmap_summary.txt`
- `out/sites_in_bins.txt`
- `out/columns/{summary.txt,spotcheck.txt,GCA_056209055.1.tsv}`
