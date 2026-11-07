"""Analyse the second kangaroo-read pass (results/trna_hunt/run_hunt2.sh).

1. CTA tRNAs: ARAGORN on every read in tloop_*.fq; keep tRNAs with a CTA (amber) anticodon; cluster identical
   sequences; for each, count reads, and how many of those reads have a mate (same spot) mapped to the kangaroo MAG
   (mag_*.sam). Mate-on-MAG = the tRNA sits within a few hundred bp of MAG sequence, i.e. in this organism's genome.
   For comparison, the same mate-linkage rate is reported for all other tRNAs found.
2. TAG-site read support: for every in-frame TAG in the MAG's genes (table 15), the codon seen in each mapped read.

usage: python src/hunt2_analyze.py
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pyaragorn
import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import read_fasta, revcomp

H = Path(__file__).resolve().parent.parent / "results" / "trna_hunt"


def spot(name):
    return name.split()[0].lstrip("@").split("/")[0]


def mag_mapped_spots():
    spots = set()
    for p in sorted(H.glob("mag_*.sam")):
        for line in open(p):
            if line[0] == "@":
                continue
            f = line.split("\t", 3)
            if len(f) > 2 and not int(f[1]) & 4:
                spots.add(f[0])
    return spots


def trnas(mapped):
    finder = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
    cta = defaultdict(list)
    n_trna = n_linked = 0
    for p in sorted(H.glob("tloop_*.fq")):
        lines = open(p).read().split("\n")
        for i in range(0, len(lines) - 3, 4):
            name, seq = lines[i], lines[i + 1]
            if not name.startswith("@"):
                continue
            for g in finder.find_rna(seq.encode()):
                n_trna += 1
                linked = spot(name) in mapped
                n_linked += linked
                if g.anticodon.upper() == "CTA":
                    cta[g.sequence().upper()].append((spot(name), linked))
    print(f"tRNAs in T-loop reads: {n_trna}; with a mate/read on the kangaroo MAG: {n_linked} "
          f"({100 * n_linked / max(n_trna, 1):.1f}%)")
    print(f"CTA-anticodon tRNA reads: {sum(len(v) for v in cta.values())} in {len(cta)} distinct sequences")
    for s, v in sorted(cta.items(), key=lambda kv: -len(kv[1]))[:25]:
        print(f"  reads {len(v):4d}  MAG-linked {sum(l for _, l in v):3d}  {s}")
    return cta


def tag_support():
    contigs = dict(read_fasta(H / "kang_mag.fa.gz")) if (H / "kang_mag.fa.gz").exists() else None
    if contigs is None:
        contigs = {}
        name = None
        for l in open(H / "kang_mag.fa"):
            if l.startswith(">"):
                name = l[1:].split()[0]
                contigs[name] = []
            else:
                contigs[name].append(l.strip())
        contigs = {k: "".join(v).upper() for k, v in contigs.items()}
    finder = pyrodigal.GeneFinder()
    finder.train(*[s.encode() for s in contigs.values()], translation_table=15)
    site_of = {}  # (contig, pos of first base in forward coords) -> (strand, site id)
    sites = []
    for c, s in contigs.items():
        for g in finder.find_genes(s.encode()):
            nt = g.sequence()
            for i in range(0, len(nt) - 3, 3):
                if nt[i:i + 3] == "TAG":
                    if g.strand == 1:
                        p0 = g.begin - 1 + i
                    else:
                        p0 = g.end - 3 - i  # forward-strand start of the CTA that is TAG on the minus strand
                    site_of[(c, p0)] = (g.strand, len(sites))
                    sites.append((c, p0, g.strand))
    obs = defaultdict(Counter)
    starts = defaultdict(list)
    for (c, p0), v in site_of.items():
        starts[c].append(p0)
    for c in starts:
        starts[c].sort()
    import bisect
    for p in sorted(H.glob("mag_*.sam")):
        for line in open(p):
            if line[0] == "@":
                continue
            f = line.split("\t")
            if len(f) < 11 or int(f[1]) & 4 or f[2] not in starts:
                continue
            c, pos, cigar, seq = f[2], int(f[3]) - 1, f[5], f[9]
            # reference -> read map for M blocks
            refmap = {}
            r, q = pos, 0
            for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                ln = int(ln)
                if op in "M=X":
                    for k in range(ln):
                        refmap[r + k] = q + k
                    r += ln
                    q += ln
                elif op in "IS":
                    q += ln
                elif op in "DN":
                    r += ln
            lo = bisect.bisect_left(starts[c], pos)
            hi = bisect.bisect_right(starts[c], r - 3)
            for p0 in starts[c][lo:hi]:
                if all(p0 + k in refmap for k in range(3)):
                    cod = "".join(seq[refmap[p0 + k]] for k in range(3))
                    strand, sid = site_of[(c, p0)]
                    if strand == -1:
                        cod = revcomp(cod)
                    obs[sid][cod] += 1
    covered = [sid for sid in range(len(sites)) if sum(obs[sid].values()) >= 3]
    fixed = sum(1 for sid in covered if obs[sid]["TAG"] / sum(obs[sid].values()) >= 0.9)
    poly = [sid for sid in covered if 0.1 < obs[sid]["TAG"] / sum(obs[sid].values()) < 0.9]
    absent = [sid for sid in covered if obs[sid]["TAG"] / sum(obs[sid].values()) <= 0.1]
    depth = sorted(sum(obs[sid].values()) for sid in covered)
    print(f"in-frame TAG sites in MAG genes: {len(sites)}; covered by >= 3 reads: {len(covered)} "
          f"(median depth {depth[len(depth) // 2] if depth else 0})")
    print(f"  TAG in >= 90% of reads: {fixed}; mixed (10-90%): {len(poly)}; TAG in <= 10%: {len(absent)}")
    alt = Counter()
    for sid in poly + absent:
        for cod, k in obs[sid].items():
            if cod != "TAG":
                alt[cod] += k
    print("  other codons at mixed/absent sites:", alt.most_common(10))


if __name__ == "__main__":
    mapped = mag_mapped_spots()
    print(f"spots with a read on the kangaroo MAG: {len(mapped)}")
    trnas(mapped)
    tag_support()
