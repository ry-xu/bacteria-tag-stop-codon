"""Re-check the claimant's per-site read tables (results/read_support/*.sites.tsv) for the candidate.
1. Coordinate/strand check: does the reference (GCA_056209055.1) read TAG at (contig, pos0, strand)?
2. Pool the four animals; list polymorphic sites (another codon in >= 2 reads of one run) and annotate each with
   (a) the residue a relative's protein has there (my tblastn sites incl. a looser alignment filter),
   (b) whether the TAG is the terminal stop of a table-11 gene (C-terminal read-through) or internal,
   (c) the codon each relative-genome ortholog alignment puts opposite it.
3. Compare: the alternative-codon spectrum at TAG sites vs the sequencing-error spectrum at TGG sites (from the
   claimant's run summaries), single-read vs >= 2-read alternatives.
usage: read_sites_check.py"""
import re
from collections import Counter, defaultdict
from pathlib import Path

D = Path(__file__).resolve().parent.parent
RS = D.parent / "results" / "read_support"
RUNS = {"CRR939056": "TA20190005", "CRR939061": "TA20190013", "CRR939060": "TA20190011", "CRR939058": "TA20190008"}
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


dna = dict(read_fasta(D / "data" / "GCA_056209055.1.fna"))
# 1. coordinate check
obs = defaultdict(dict)
strand_of = {}
bad = Counter()
for run in RUNS:
    for line in list(open(RS / f"GCA_056209055.1_{run}.sites.tsv"))[1:]:
        c, p, st, cods = line.rstrip("\n").split("\t")
        p, st = int(p), int(st)
        ref = dna[c][p:p + 3]
        ref = rc(ref) if st == -1 else ref
        bad[ref == "TAG"] += 1
        strand_of[(c, p)] = st
        obs[(c, p)][run] = Counter({k: int(v) for k, v in (x.rsplit(":", 1) for x in cods.split())})
print(f"coordinate/strand check: reference codon == TAG at {bad[True]} site-rows, != TAG at {bad[False]}")
print(f"  strands of sites: {Counter(strand_of.values())}")

# gene ends (table 11 pyrodigal) for the candidate
gene_end = {}
for line in list(open(D / "out" / "genes" / "GCA_056209055.1.tsv"))[1:]:
    gid, c, b, e, st, stop, part = line.rstrip("\n").split("\t")
    b, e, st = int(b), int(e), int(st)
    p0 = e - 3 if st == 1 else b - 1
    gene_end[(c, p0)] = gid
# my tblastn TAG sites (strict and loose)
strict = {}
for f in ["cand_all.sites.tsv"]:
    p = D / "out" / "stops" / f
    if p.exists():
        for line in list(open(p))[1:]:
            x = line.rstrip("\n").split("\t")
            strict[(x[0], int(x[1]))] = (x[3], x[4], x[5].split("|")[0][:15], x[6], x[7])
loose = {}
p = D / "out" / "stops" / "cand_loose.sites.tsv"
if p.exists():
    for line in list(open(p))[1:]:
        x = line.rstrip("\n").split("\t")
        loose[(x[0], int(x[1]))] = (x[3], x[4], x[5].split("|")[0][:15], x[6], x[7])

# 2. pooled polymorphic sites
rows = []
for k, per in obs.items():
    alts = Counter()
    for run, o in per.items():
        for cod, n in o.items():
            if cod != "TAG" and n >= 2:
                alts[cod] += 1
    if alts:
        rows.append((k, per, alts))
print(f"\npolymorphic sites (alt codon in >= 2 reads of one run): {len(rows)}; alt codons per site: "
      f"{Counter(a for _, _, al in rows for a in al)}")
print("contig\tpos0\tstrand\tper-run codons\tgene_end(table11)\trelative residue (strict | loose)")
for k, per, alts in sorted(rows):
    st = strand_of[k]
    pr = " | ".join(f"{RUNS[r][-2:]}:" + ",".join(f"{c}{n}" for c, n in o.most_common()) for r, o in per.items())
    ge = "END" if k in gene_end else "internal"
    print(f"{k[0]}\t{k[1]}\t{st}\t{pr}\t{ge}\t{strict.get(k, '-')} | {loose.get(k, '-')}")

# TAG sites overall: how many are table-11 gene ends vs internal; how many have a relative residue
tot = Counter(("END" if k in gene_end else "internal") for k in obs)
print(f"\nall TAG sites with reads: {len(obs)}; table-11 gene terminal stops: {tot['END']}, internal: {tot['internal']}")
print(f"  with a relative's residue (strict tblastn filter): {sum(k in strict for k in obs)}; loose: {sum(k in loose for k in obs)}")
fixed_all = [k for k, per in obs.items() if all(o['TAG'] / sum(o.values()) >= 0.9 for o in per.values() if sum(o.values()) >= 3) and any(sum(o.values()) >= 3 for o in per.values())]
print(f"  fixed in every covering run: {len(fixed_all)}; of those END {sum(k in gene_end for k in fixed_all)}")
poly = {k for k, _, _ in rows}
print(f"  polymorphic: END {sum(k in gene_end for k in poly)} internal {sum(k not in gene_end for k in poly)}")
# residue at polymorphic vs fixed (loose)
for name, ks in (("polymorphic", poly), ("non-polymorphic", set(obs) - poly)):
    print(f"  {name}: relative residue (loose) {Counter(loose[k][1] for k in ks if k in loose).most_common()}")

# 3. single-read vs multi-read alternatives at TAG sites
one, two = Counter(), Counter()
for k, per in obs.items():
    for run, o in per.items():
        for cod, n in o.items():
            if cod == "TAG":
                continue
            (two if n >= 2 else one)[cod] += 1
print(f"\nTAG sites, alt codon seen in exactly 1 read of a run (site-runs): {one.most_common()}")
print(f"TAG sites, alt codon seen in >= 2 reads of a run (site-runs): {two.most_common()}")
