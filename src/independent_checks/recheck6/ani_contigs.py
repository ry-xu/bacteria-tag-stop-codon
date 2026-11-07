"""Per-contig and genome-level identity of a query bin vs a reference (minimap2 -c PAF, primary only).
usage: ani_contigs.py <query.fna> <paf> <label>  -> prints summary; writes out/<label>.contigs.tsv"""
import sys
from collections import defaultdict
from pathlib import Path

fa, paf, label = sys.argv[1:4]
lens, name = {}, None
for line in open(fa):
    if line[0] == ">":
        name = line[1:].split()[0]; lens[name] = 0
    else:
        lens[name] += len(line.strip())
cov = defaultdict(list); mat = defaultdict(int); blk = defaultdict(int)
for line in open(paf):
    f = line.split("\t")
    if "tp:A:P" not in line:
        continue
    q, qs, qe = f[0], int(f[2]), int(f[3])
    nm = int(next(t[5:] for t in f[12:] if t.startswith("NM:i:")))
    cov[q].append((qs, qe)); mat[q] += int(f[9]); blk[q] += int(f[10])
def covered(iv):
    iv = sorted(iv); tot = 0; cs = ce = -1
    for s, e in iv:
        if s > ce:
            tot += ce - cs if ce > cs else 0; cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if ce > cs else 0)
out = Path(__file__).resolve().parent.parent / "out" / f"{label}.contigs.tsv"
rows = []
for c, L in lens.items():
    cv = covered(cov[c]) if cov[c] else 0
    rows.append((c, L, cv, mat[c] / blk[c] if blk[c] else 0))
with open(out, "w") as fh:
    fh.write("contig\tlen\taligned\tident\n")
    for r in rows:
        fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\n")
G = sum(lens.values()); A = sum(r[2] for r in rows)
ident = sum(mat.values()) / max(1, sum(blk.values()))
bins = {"none(<10%)": 0, "partial": 0, ">=80% aligned": 0}
for r in rows:
    fr = r[2] / r[1]
    bins["none(<10%)" if fr < 0.1 else (">=80% aligned" if fr >= 0.8 else "partial")] += r[1]
hi = sum(r[1] for r in rows if r[2] / r[1] >= 0.5 and r[3] >= 0.95)
print(f"{label}: {G} bp, {len(lens)} contigs; aligned {A} bp ({A/G:.1%}); identity in aligned blocks {ident:.4f}; "
      f"bp by contig class: {bins}; bp in contigs >=50% aligned at >=95% id: {hi} ({hi/G:.1%})")
