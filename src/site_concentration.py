"""Are read-through sites of a protein family concentrated at one position more than chance allows? For each
(stop codon, protein family) the read-through hits (one per genome, gap-free windows only, genome-wide readers left
out as in site_clusters.py) are binned by subject position (15 residues); under the null they fall uniformly over the
bins the family's hits span. Score = binomial upper tail of the fullest bin x number of bins (Bonferroni within the
family). Sec/Pyl sites give astronomically small values; scattered pseudogene stops don't.
Fusion junctions (two adjacent genes whose Swiss-Prot homolog is one fused protein) concentrate too, but at every stop
codon in proportion to its use; recoding is stop-specific. stop_share = this stop's hits / all stops' hits at the bin;
"specific" = stop_share >= 0.8. Split genes whose stop codon is inherited within one clade are stop-specific too, so\nalso require the top bin to span >= 3 GTDB orders (oversampled clades: Aristaeellaceae, Collinsella).
usage: python src/site_concentration.py  -> results/site_concentration.tsv (sorted by p)"""
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

from site_clusters import BIN, family

ROOT = Path(__file__).resolve().parent.parent


def logbinom_tail(k, n, p):
    """log10 P(X >= k), X ~ Bin(n, p)."""
    terms = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(p) + (n - i) * math.log1p(-p)
             for i in range(k, n + 1)]
    m = max(terms)
    return (m + math.log(sum(math.exp(t - m) for t in terms))) / math.log(10)


def main():
    wide = set()
    for f in (ROOT / "results" / "rt").glob("*.stats.json"):
        d = json.loads(f.read_text())
        rt, te, rates = d["readthrough"], d["tested"], d.get("rates") or {}
        est = {x: rt[x] / rates.get(x, 1.0) for x in rt}
        for s in rt:
            others = max([est[o] for o in est if o != s] + [1])
            if est[s] >= 20 and est[s] / others >= 5 and rt[s] / max(1, te[s]) >= 0.1:
                wide.add((d["acc"], s))
    tax = {}
    for dom in ("ar53", "bac120"):
        with gzip.open(ROOT / "data" / f"{dom}_metadata.tsv.gz", "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                tax[r["accession"][3:]] = r["gtdb_taxonomy"].split(";")
    order = lambda a: tax[a][3] if a in tax else "?"  # genomes newer than GTDB (mostly Lentisphaeria) count as one lineage
    hits = defaultdict(dict)  # (stop, family) -> {acc: (bin, sres, sseqid)}
    for f in (ROOT / "results" / "rt").glob("*.tsv"):
        for r in csv.DictReader(open(f), delimiter="\t"):
            if (r["acc"], r["stop"]) in wide or r["spos"] in ("None", "") or int(r.get("gaps_win") or 1) > 0:
                continue
            hits[(r["stop"], family(r["stitle"]))][r["acc"]] = (int(r["spos"]) // BIN, r["sres"], r["sseqid"])
    rows = []
    for (stop, fam), g in hits.items():
        n = len(g)
        if n < 5:
            continue
        bins, accs = defaultdict(list), defaultdict(list)
        for acc, (b, sres, sid) in g.items():
            bins[b].append(sres)
            accs[b].append(acc)
        nb = max(bins) - min(bins) + 1
        top = max(bins, key=lambda b: len(bins[b]))
        k = len(bins[top])
        if nb < 2:
            continue
        lp = logbinom_tail(k, n, 1 / nb) + math.log10(nb)
        sres = max(set(bins[top]), key=bins[top].count)
        at_bin = {x: sum(1 for v in hits.get((x, fam), {}).values() if v[0] == top) for x in ("TAA", "TAG", "TGA")}
        rows.append(dict(stop=stop, family=fam, genomes=n, top_bin=f"{top * BIN}-{top * BIN + BIN - 1}", in_top=k,
                         TAA_TAG_TGA_at_bin="/".join(str(at_bin[x]) for x in ("TAA", "TAG", "TGA")),
                         stop_share=round(k / max(1, sum(at_bin.values())), 2),
                         bins_spanned=nb, log10p=round(lp, 1), subject_residue=sres,
                         orders=len({order(a) for a in accs[top]} - {"?"}),
                         phyla=len({tax[a][1] if a in tax else "?" for a in accs[top]}),
                         share=round(bins[top].count(sres) / k, 2)))
    rows.sort(key=lambda x: x["log10p"])
    cols = list(rows[0])
    with open(ROOT / "results" / "site_concentration.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for x in rows:
            f.write("\t".join(str(x[c]) for c in cols) + "\n")
    m = len(rows)
    print(f"{m} (stop, family) pairs with >= 5 genomes; Bonferroni over pairs: log10 p < {math.log10(0.01 / m):.1f}")
    for stop in ("TAG", "TAA", "TGA"):
        sig = [x for x in rows if x["stop"] == stop and x["log10p"] < math.log10(0.01 / m)]
        spec = [x for x in sig if x["stop_share"] >= 0.8 and x["orders"] >= 3]
        print(f"\n== {stop}: {len(sig)} concentrated families, {len(spec)} of them specific to {stop} and in >= 3 orders")
        sig = spec
        for x in sig[:40]:
            print(f"  {x['family'][:55]:55s} {x['top_bin']:>9} {x['in_top']:4d}/{x['genomes']:<5d} bins {x['bins_spanned']:4d} "
                  f"log10p {x['log10p']:7.1f} ord {x['orders']:3d} phy {x['phyla']:2d} stops {x['TAA_TAG_TGA_at_bin']:>10} subj {x['subject_residue']}({x['share']})")


if __name__ == "__main__":
    main()
