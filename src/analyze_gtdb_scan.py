"""Summarise results/gtdb_stops_rf.tsv (stream scan) with GTDB taxonomy and CheckM2 quality.

Candidates for a reassigned stop codon:
  noRF1  bacteria with no RF1 hit (RF1 reads UAG), completeness >= MINC   -> UAG free to be reassigned
  noRF2  bacteria with no RF2 hit (RF2 reads UGA)                          -> UGA (known: code 4 Mycoplasmatales etc.)
  rareTAG / rareTGA / rareTAA: a stop codon ending < 0.5% of complete genes
A missing gene in a MAG can be incompleteness, so each candidate gets P(absent by chance) ~ (1 - completeness), and
lineages are summarised at genus level (several independent genomes lacking RF1 is much stronger than one).

usage: python src/analyze_gtdb_scan.py [min_completeness=90]
"""
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MINC = float(sys.argv[1]) if len(sys.argv) > 1 else 90.0


def meta():
    m = {}
    for dom in ("ar53", "bac120"):
        with gzip.open(ROOT / "data" / f"{dom}_metadata.tsv.gz", "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["gtdb_representative"] == "t":
                    m[r["accession"]] = dict(tax=r["gtdb_taxonomy"], comp=float(r["checkm2_completeness"] or 0),
                                             cont=float(r["checkm2_contamination"] or 0), gc=r["gc_percentage"],
                                             cat=r["ncbi_genome_category"])
    return m


def main():
    M = meta()
    rows = list(csv.DictReader(open(ROOT / "results" / "gtdb_stops_rf.tsv"), delimiter="\t"))
    bac = [r for r in rows if r["domain"] == "bacteria"]
    print(f"scanned {len(rows)} genomes ({len(bac)} bacteria)")
    out = {k: [] for k in ("noRF1", "noRF2", "rareTAG", "rareTGA", "rareTAA")}
    for r in bac:
        q = M.get(r["acc"])
        if not q:
            continue
        stops = {k: int(r[k]) for k in ("TAA", "TAG", "TGA")}
        n = sum(stops.values()) or 1
        rec = dict(acc=r["acc"], tax=q["tax"], comp=q["comp"], cont=q["cont"], gc=q["gc"], **stops,
                   fTAG=stops["TAG"] / n, fTGA=stops["TGA"] / n, fTAA=stops["TAA"] / n, RF1=r["n_RF1"], RF2=r["n_RF2"],
                   markers=r.get("markers", ""))
        if q["comp"] < MINC:
            continue
        if r["n_RF1"] == "0":
            out["noRF1"].append(rec)
        if r["n_RF2"] == "0":
            out["noRF2"].append(rec)
        for k in ("TAG", "TGA", "TAA"):
            if rec["f" + k] < 0.005:
                out["rare" + k].append(rec)
    hq = sum(1 for r in bac if M.get(r["acc"], {}).get("comp", 0) >= MINC)
    print(f"bacteria with completeness >= {MINC}: {hq}")
    for k, recs in out.items():
        print(f"\n== {k}: {len(recs)} genomes")
        by = defaultdict(list)
        for x in recs:
            by[";".join(x["tax"].split(";")[1:6])].append(x)
        for g, xs in sorted(by.items(), key=lambda kv: -len(kv[1]))[:40]:
            ex = xs[0]
            print(f"  {len(xs):3d}  {g.replace('p__','').replace('c__','').replace('o__','').replace('f__','').replace('g__','')}"
                  f"  e.g. {ex['acc']} comp={ex['comp']:.0f} TAA/TAG/TGA={ex['TAA']}/{ex['TAG']}/{ex['TGA']} RF1={ex['RF1']} RF2={ex['RF2']}")
        with open(ROOT / "results" / f"gtdb_{k}.tsv", "w") as f:
            cols = list(recs[0].keys()) if recs else []
            f.write("\t".join(cols) + "\n")
            for x in recs:
                f.write("\t".join(str(x[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
