"""Anticodon-stem pairing of the tRNAs that are near-cognate to UAG, readers vs relatives.

Idea (Kachale et al. 2023, Nature; Fakih et al. 2026, NAR): a tRNA whose anticodon stem has only 4 paired bases
instead of 5 can read a stop codon that differs from its own codon at one position (Blastocrithidia's tRNA-Trp(CCA)
reads UGA). The Victivallales TAG readers have no CUA-anticodon tRNA, so if one of their tRNAs one mismatch away from
UAG (Lys CUU, Gln CUG, Glu CUC, Ser CGA, Trp CCA, Leu CAA, Tyr GUA, or their U34 versions) has a 4-bp stem that its
relatives' copies lack, that tRNA is the candidate decoder and its amino acid the candidate meaning.

Stem = the 5 base pairs flanking the 7-nt anticodon loop (positions 27-31 : 39-43); Watson-Crick and G-U count.
usage: python src/trna_stems.py <acc_list: acc<TAB>label> -> results/trna_stems.tsv + summary
"""
import sys
from collections import defaultdict
from pathlib import Path

import pyaragorn

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import ROOT, fetch, read_fasta

PAIRS = {("G", "C"), ("C", "G"), ("A", "T"), ("T", "A"), ("G", "T"), ("T", "G")}
NEAR = {"CTT", "TTT", "CTG", "TTG", "CTC", "TTC", "CGA", "TGA", "CCA", "CAA", "TAA", "GTA", "ATA", "CTA"}


def stem(seq, off):
    a = off - 1  # 0-based start of the anticodon
    if a - 7 < 0 or a + 9 >= len(seq):
        return None
    pairs = [(seq[a - 7 + i], seq[a + 9 - i]) for i in range(5)]  # 27:43 ... 31:39
    ok = [p in PAIRS for p in pairs]
    return "".join("|" if x else "." for x in ok)  # top (27:43) first


def main(listfile):
    rows = []
    for line in open(listfile):
        if not line.strip():
            continue
        acc, label = line.split()[:2]
        finder = pyaragorn.RNAFinder(11, tmrna=False)
        for c, s in read_fasta(fetch(acc)):
            for g in finder.find_rna(s.encode()):
                ac = g.anticodon.upper().replace("U", "T")
                if g.type != "tRNA" or ac not in NEAR:
                    continue
                seq = g.sequence().upper()
                st = stem(seq, g.anticodon_offset)
                rows.append((label, acc, g.amino_acid, ac, st, st.count("|") if st else -1, seq))
    out = ROOT / "results" / "trna_stems.tsv"
    with open(out, "w") as f:
        f.write("label\tacc\taa\tanticodon\tstem\tpairs\tseq\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    # summary: per label x (aa, anticodon), how many genomes have a copy with < 5 stem pairs, of genomes with any copy
    have, short = defaultdict(set), defaultdict(set)
    for label, acc, aa, ac, st, n, seq in rows:
        have[(label, aa, ac)].add(acc)
        if 0 <= n < 5:
            short[(label, aa, ac)].add(acc)
    labels = sorted({r[0] for r in rows})
    keys = sorted({(aa, ac) for _, aa, ac in have})
    print("genomes with a <5-bp anticodon stem copy / genomes with the tRNA")
    print("tRNA".ljust(10) + "".join(l.ljust(16) for l in labels))
    for aa, ac in keys:
        print(f"{aa}-{ac}".ljust(10) + "".join(
            f"{len(short[(l, aa, ac)])}/{len(have[(l, aa, ac)])}".ljust(16) for l in labels))


if __name__ == "__main__":
    main(sys.argv[1])
