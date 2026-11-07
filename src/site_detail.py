"""Taxonomic spread and gene context of one read-through site cluster (stop, family, subject-position bin).
usage: python src/site_detail.py <stop> "<family substring>" <bin_start>"""
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from site_clusters import BIN, family

ROOT = Path(__file__).resolve().parent.parent
stop, fam_sub, b0 = sys.argv[1], sys.argv[2].lower(), int(sys.argv[3])
tax = {}
for dom in ("ar53", "bac120"):
    with gzip.open(ROOT / "data" / f"{dom}_metadata.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tax[r["accession"][3:]] = r["gtdb_taxonomy"]
rows = []
for f in (ROOT / "results" / "rt").glob("*.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["stop"] != stop or r["spos"] in ("None", "") or int(r.get("gaps_win") or 1) > 0:
            continue
        if fam_sub in family(r["stitle"]) and int(r["spos"]) // BIN * BIN == b0:
            rows.append(r)
print(len(rows), "hits")
lv = lambda a, i: (tax.get(a, "?").split(";") + ["?"] * 7)[i]
for i, name in ((1, "phylum"), (3, "order"), (4, "family")):
    c = Counter(lv(r["acc"], i) for r in rows)
    print(f"{name}: {len(c)} ->", c.most_common(8))
print("subject:", Counter(r["sseqid"].split("|")[-1] for r in rows).most_common(5))
print("subject residue at stop:", Counter(r["sres"] for r in rows).most_common(5))
print("next in-frame gene starts (codons after stop; -1 none):", Counter(int(r["next_gene_codons"]) for r in rows).most_common(6))
print("ext reaches contig end:", Counter(r["ext_to_contig_end"] for r in rows))
print("query len - up_len (aa after stop, incl. X):", sorted(int(r["len"]) - int(r["up_len"]) for r in rows)[::max(1, len(rows) // 10)])
print("pident:", sorted(float(r["pident"]) for r in rows)[::max(1, len(rows) // 10)])
for r in rows[:4]:
    print(" ", r["acc"], lv(r["acc"], 4), r["gene"], r["qwin"], r["swin"])
