"""Cross-genome clustering of read-through sites: the Sec/Pyl signature is the same in-frame stop at the same site of
the same protein family in many genomes; pseudogenes and sequencing errors are scattered.

Input: results/rt/*.tsv (one row per read-through gene: stop codon, best Swiss-Prot subject, subject position
aligned to the stop, subject residue there). A "site" = (stop codon, protein family, subject position // BIN).
Family = Swiss-Prot protein name with strain/number suffixes stripped, so homologous subjects pool together;
positions from different subjects of one family are only approximately comparable, hence the coarse bin.
Genomes with genome-wide read-through at a stop (>= 20 sampling-corrected read-throughs, >= 5x the other stops,
rate >= 0.1) are left out for that stop.
Only hits whose alignment is gap-free within 12 columns of the stop are used (new-format results only): gene pairs
whose Swiss-Prot homolog is a natural fusion (feoA/feoB, lpxC/fabZ, pgk/tpi ...) otherwise dominate.

usage: python src/site_clusters.py [min_genomes=3]  -> results/site_clusters.tsv
"""
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = 15
MING = int(sys.argv[1]) if __name__ == "__main__" and len(sys.argv) > 1 else 3


def taxonomy():
    t = {}
    for dom in ("ar53", "bac120"):
        with gzip.open(ROOT / "data" / f"{dom}_metadata.tsv.gz", "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                t[r["accession"][3:]] = r["gtdb_taxonomy"]
    return t


def family(stitle):
    name = stitle.split(" ", 1)[1] if " " in stitle else stitle
    name = re.sub(r"\b(homolog|putative|probable|uncharacterized)\b", "", name, flags=re.I)
    name = re.sub(r"\s+\d+$|\s+[A-Z]$|\s+subunit\s+\w+$", "", name.strip())
    name = re.sub(r"(?<=[A-Za-z])\d+$", "", name.strip())  # MtbB1/MtbB3 -> MtbB, FdhA2 -> FdhA
    return name.strip().lower()


def main():
    tax = taxonomy()
    wide = set()
    for f in (ROOT / "results" / "rt").glob("*.stats.json"):
        d = json.loads(f.read_text())
        rt, te, rates = d["readthrough"], d["tested"], d.get("rates") or {}
        est = {x: rt[x] / rates.get(x, 1.0) for x in rt}  # subsampled runs: scale counts back up
        for s in ("TAG", "TGA", "TAA"):
            others = max([est[o] for o in est if o != s] + [1])
            if est[s] >= 20 and est[s] / others >= 5 and rt[s] / max(1, te[s]) >= 0.1:
                wide.add((d["acc"], s))
    print("genome-wide read-through (excluded from site clusters):", len(wide), "genome x stop")
    sites = defaultdict(list)
    for f in (ROOT / "results" / "rt").glob("*.tsv"):
        for r in csv.DictReader(open(f), delimiter="\t"):
            if (r["acc"], r["stop"]) in wide or r["spos"] in ("None", ""):
                continue
            # fusion junctions and frameshift artifacts leave gaps next to the stop; Sec/Pyl sites align without gaps
            if "gaps_win" not in r or int(r["gaps_win"]) > 0:
                continue
            key = (r["stop"], family(r["stitle"]), int(r["spos"]) // BIN)
            sites[key].append(r)
    out = []
    for (stop, fam, b), rs in sites.items():
        genomes = {r["acc"] for r in rs}
        if len(genomes) < MING:
            continue
        lineages = {";".join(tax.get(a, "?").split(";")[:6]) for a in genomes}
        fams = {";".join(tax.get(a, "?").split(";")[:5]) for a in genomes}
        sres = defaultdict(int)
        for r in rs:
            sres[r["sres"]] += 1
        top = max(sres, key=sres.get)
        out.append(dict(stop=stop, family=fam, pos=f"{b * BIN}-{b * BIN + BIN - 1}", genomes=len(genomes),
                        genera=len(lineages), families=len(fams), subject_residue=top,
                        subject_residue_share=round(sres[top] / len(rs), 2),
                        example=rs[0]["sseqid"], phyla=";".join(sorted({(tax.get(a, "?;?").split(";") + ["?"])[1] for a in genomes}))))
    # flags: known = Sec/Pyl in the subject (U/O show as X), or TGA aligned to Cys (Sec/Cys homolog pairs);
    # multistop = the same family+position also clusters under another stop codon (fusion/frameshift artifact)
    stops_at = defaultdict(set)
    for x in out:
        stops_at[(x["family"], x["pos"])].add(x["stop"])
    for x in out:
        if x["subject_residue"] == "X" or (x["stop"] == "TGA" and x["subject_residue"] == "C"):
            x["flag"] = "known"
        elif len(stops_at[(x["family"], x["pos"])]) > 1:
            x["flag"] = "multistop"
        else:
            x["flag"] = "candidate"
    known_fams = {(x["family"], x["stop"]) for x in out if x["flag"] == "known"}
    for x in out:
        if x["flag"] == "candidate" and (x["family"], x["stop"]) in known_fams:
            x["flag"] = "knownfam"  # a known Sec/Pyl family whose site fell in a neighbouring bin
    out.sort(key=lambda x: (x["flag"] != "candidate", -x["families"], -x["genera"], -x["genomes"]))
    cols = ["flag", "stop", "family", "pos", "genomes", "genera", "families", "subject_residue",
            "subject_residue_share", "example", "phyla"]
    with open(ROOT / "results" / "site_clusters.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for x in out:
            f.write("\t".join(str(x[c]) for c in cols) + "\n")
    print(f"{len(out)} site clusters with >= {MING} genomes")
    from collections import Counter as _C
    print("flags:", dict(_C(x["flag"] for x in out)))
    for x in out[:30]:
        print(f"  {x['flag'][:5]} {x['stop']} {x['family'][:45]:45s} {x['pos']:>9} genomes={x['genomes']:3d} genera={x['genera']:3d} "
              f"fams={x['families']:3d} subj={x['subject_residue']}({x['subject_residue_share']}) {x['phyla'][:50]}")


if __name__ == "__main__":
    main()
