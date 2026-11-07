"""Pull every tRNA with a stop anticodon (CTA/TCA/TTA, including ARAGORN 'SeC') out of codon-hunt's archived
genome-wide ARAGORN screen (191,958 GTDB R232 species reps; read-only archive) and attach GTDB taxonomy.

usage: TRNA_SCREEN_DIR=<dir> uv run python src/extract_stop_trnas.py  ->  results/stop_trnas.tsv, results/genomes_screened.tsv
"""
import csv
import gzip
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# NOT PORTABLE: needs the archived ARAGORN screen of all GTDB R232 species reps (trnas_*.tsv.gz), which is not
# public. Point TRNA_SCREEN_DIR at an equivalent screen to rerun; results/stop_trnas.tsv is not shipped either.
ARCH = Path(os.environ.get("TRNA_SCREEN_DIR", ROOT / "data" / "trna_screen"))
STOP = {"CTA": "TAG", "TCA": "TGA", "TTA": "TAA"}


def taxonomy():
    tax = {}
    for dom in ("ar53", "bac120"):
        with gzip.open(ROOT / "data" / f"{dom}_metadata.tsv.gz", "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["gtdb_representative"] == "t":
                    tax[r["accession"][3:]] = (r["gtdb_taxonomy"], r["checkm2_completeness"], r["checkm2_contamination"],
                                               r["gc_percentage"], r["genome_size"], r["ncbi_genome_category"])
    return tax


def main():
    tax = taxonomy()
    out = open(ROOT / "results" / "stop_trnas.tsv", "w")
    w = None
    n = 0
    for p in sorted(ARCH.glob("trnas_*.tsv.gz")):
        with gzip.open(p, "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["anticodon"] in STOP:
                    t = tax.get(r["acc"], ("",) * 6)
                    r.update(codon=STOP[r["anticodon"]], taxonomy=t[0], completeness=t[1], contamination=t[2],
                             gc=t[3], size=t[4], category=t[5])
                    if w is None:
                        w = csv.DictWriter(out, fieldnames=list(r), delimiter="\t")
                        w.writeheader()
                    w.writerow(r)
                    n += 1
    out.close()
    print("stop-anticodon tRNAs:", n)


if __name__ == "__main__":
    main()
