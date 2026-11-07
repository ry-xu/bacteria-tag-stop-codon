"""Homology-free check: for table-11 genes ending in TAG vs TAA/TGA, how far does the reading frame stay open after
the stop (codons to the next in-frame TAA or TGA, TAGs ignored)? A genome that reads TAG as sense has many
'genes' truncated at a TAG with a long open frame behind it; a TAG-stopping genome does not.
usage: downstream_orf.py <genome> ..."""
import sys
from collections import defaultdict
from pathlib import Path

D = Path(__file__).resolve().parent.parent
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(p):
    name, seq = None, []
    for line in open(p):
        if line[0] == ">":
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip().upper())
    if name:
        yield name, "".join(seq)


print("genome\tstop\tn_genes\tfrac_open>=50codons\tfrac_open>=100codons\tmedian_open")
for g in sys.argv[1:]:
    dna = dict(read_fasta(D / "data" / f"{g}.fna"))
    dist = defaultdict(list)
    for line in list(open(D / "out" / "genes" / f"{g}.tsv"))[1:]:
        gid, c, b, e, st, stop, part = line.rstrip("\n").split("\t")
        if part != "00" or stop == "NA":
            continue
        b, e, st = int(b), int(e), int(st)
        s = dna[c]
        down = s[e:] if st == 1 else rc(s[:b - 1])
        n, closed = 0, False
        for k in range(0, len(down) - 2, 3):
            cod = down[k:k + 3]
            if cod in ("TAA", "TGA"):
                closed = True
                break
            n += 1
        if not closed:
            continue          # runs off the contig: unknown
        dist[stop].append(n)
    for stop in ("TAG", "TAA", "TGA"):
        v = sorted(dist[stop])
        if not v:
            continue
        print(f"{g}\t{stop}\t{len(v)}\t{sum(x >= 50 for x in v) / len(v):.2f}\t{sum(x >= 100 for x in v) / len(v):.2f}\t{v[len(v) // 2]}")
