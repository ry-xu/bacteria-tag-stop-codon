"""Rank every genome tested so far (results/rt/*.stats.json, any sampling rates) for genome-wide read-through of one
stop codon. Subsampled screens test different fractions of each stop, so compare rates, not counts:
  rate(stop) = genes read through at that stop / genes tested at that stop
  expected   = tested(stop) * pooled rate at the other two stops (+ 0.5 pseudo-count)
  excess     = observed / expected
A genome-wide reader (the four known TAG readers) has tens to hundreds of read-throughs at one stop and ~0 at the
others. Flag: observed >= MINOBS, excess >= MINX and rate >= MINRATE (the known readers: 25-50% of genes ending in the
reassigned stop; code-34 archaea 15-40%; gene-specific Sec/Pyl and noise: 0.5-3%).
Rates alone favour genomes that rarely end genes in the test stop: a noisy MAG with ~20 scattered read-throughs at each
stop shows 25% at a rare TAG and 2% at a common TAA (UBA11524, 2026-11-05). So also require absolute specificity:
  abs = (observed / sampling rate) / max over the other stops of (observed / sampling rate), >= MINABS (4)
(known readers: 9-56).

usage: python src/rank_readers.py [MINOBS=10] [MINX=8] [MINRATE=0.1] [MINABS=4]   -> results/readers_<stop>.tsv, top of each list printed
Reads results/rt/*.stats.json if present (your own screen), otherwise the shipped data/screen/screened_genomes.tsv.gz.
"""
import csv
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MINOBS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
MINX = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
MINRATE = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
MINABS = float(sys.argv[4]) if len(sys.argv) > 4 else 4.0


def load_stats():
    """Per-genome read-through stats: results/rt/*.stats.json if present, else the shipped screen table."""
    rt = sorted((ROOT / "results" / "rt").glob("*.stats.json"))
    if rt:
        return [json.loads(f.read_text()) for f in rt]
    out = []
    with gzip.open(ROOT / "data" / "screen" / "screened_genomes.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            st = ("TAA", "TAG", "TGA")
            rates = {s: float(r[f"sampling_{s}"]) for s in st}
            out.append(dict(acc=r["acc"], contigs=int(r["contigs"]), bp=int(r["bp"]),
                            tested={s: int(r[f"tested_{s}"]) for s in st},
                            readthrough={s: int(r[f"readthrough_{s}"]) for s in st},
                            rates=rates if any(v != 1 for v in rates.values()) else {}))
    return out


def taxonomy():
    t = {}
    for dom in ("ar53", "bac120"):
        fn = ROOT / "data" / f"{dom}_metadata.tsv.gz"
        if not fn.exists():  # GTDB metadata (scripts/setup_reference_data.sh); only used for labels
            print(f"note: {fn.name} not found, taxonomy columns left empty")
            continue
        with gzip.open(fn, "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                a = r["accession"][3:]
                t[a] = (r["gtdb_taxonomy"], r["checkm2_completeness"])
                g = r["ncbi_genbank_assembly_accession"]
                if g and g not in t:
                    t[g] = t[a]
    return t


def main():
    tax = taxonomy()
    rows = []
    for d in load_stats():
        rt, te = d["readthrough"], d["tested"]
        for s in ("TAG", "TGA", "TAA"):
            o = [x for x in ("TAA", "TAG", "TGA") if x != s]
            other_rate = sum(rt[x] for x in o) / max(1, sum(te[x] for x in o))
            exp = te[s] * other_rate + 0.5
            t, comp = tax.get(d["acc"], ("", ""))
            est = {x: rt[x] / (d.get("rates") or {}).get(x, 1.0) for x in rt}
            ab = round(est[s] / max(1.0, max(est[x] for x in o)), 1)
            rows.append(dict(stop=s, acc=d["acc"], obs=rt[s], tested=te[s], rate=round(rt[s] / max(1, te[s]), 3),
                             exp=round(exp, 2),
                             excess=round(rt[s] / exp, 1), abs=ab, other=f"{rt[o[0]]}/{te[o[0]]} {rt[o[1]]}/{te[o[1]]}",
                             sampled="y" if d.get("rates") else "n", comp=comp, tax=t))
    n = len(rows) // 3
    print(f"{n} genomes tested")
    cols = ["stop", "acc", "obs", "tested", "rate", "exp", "excess", "abs", "other", "sampled", "comp", "tax"]
    for s in ("TAG", "TGA", "TAA"):
        rs = sorted((r for r in rows if r["stop"] == s), key=lambda r: -r["excess"])
        hits = [r for r in rs if r["obs"] >= MINOBS and r["excess"] >= MINX and r["rate"] >= MINRATE and r["abs"] >= MINABS]
        (ROOT / "results").mkdir(exist_ok=True)
        with open(ROOT / "results" / f"readers_{s}.tsv", "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rs[:500]:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\n== {s}: {len(hits)} genomes with >= {MINOBS} read-throughs, >= {MINX}x the other stops' rate, "
              f"rate >= {MINRATE}, abs >= {MINABS}")
        for r in hits[:60]:
            short = ";".join(x[3:] for x in r["tax"].split(";")[1:6])
            print(f"  {r['acc']}  {r['obs']:4d}/{r['tested']:<5d} {r['rate']:<6} x{r['excess']:<6} abs {r['abs']:<5} others {r['other']:<14} "
                  f"comp {r['comp']:>5}  {short}")


if __name__ == "__main__":
    main()
