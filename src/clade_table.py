"""Table of read-through + RF status for a list of genomes with a taxonomy label.
usage: python src/clade_table.py <tsv: acc<TAB>label> """
import csv, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
rf = {r["acc"]: r for r in csv.DictReader(open(ROOT / "results" / "rf_profile.tsv"), delimiter="\t")}
rows = []
for line in open(sys.argv[1]):
    acc, label = line.rstrip("\n").split("\t")[:2]
    f = ROOT / "results" / "rt" / f"{acc}.stats.json"
    if not f.exists():
        rows.append((label, acc, "", "", "", "", "", ""))
        continue
    d = json.loads(f.read_text())
    r = rf.get(acc, {})
    rows.append((label, acc, d["readthrough"]["TAG"], d["tested"]["TAG"], d["readthrough"]["TAA"] + d["readthrough"]["TGA"],
                 r.get("frac_TAG", ""), r.get("RF1", ""), r.get("RF2", "")))
print("label\tacc\tTAG_readthrough\tTAG_tested\tTAA+TGA_readthrough\tfrac_TAG_stop\tRF1\tRF2")
for x in sorted(rows):
    print("\t".join(map(str, x)))
