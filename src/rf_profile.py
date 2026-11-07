"""Release factors and stop-codon usage per genome.
RF1 (prfA, reads UAA/UAG, PxT motif) vs RF2 (prfB, reads UAA/UGA, SPF motif), from pyrodigal proteins (table 11)
aligned to Swiss-Prot RF1/RF2; a hit is classified by its best subject. Stop usage = last codon of complete genes.
usage: python src/rf_profile.py <acc_list> [workers]  -> results/rf_profile.tsv (appends)"""
import re
import subprocess
import sys
import tempfile
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import ROOT, fetch, read_fasta

RFDB = ROOT / "data" / "rf.dmnd"
OUT = ROOT / "results" / "rf_profile.tsv"
COLS = ["acc", "bp", "contigs", "genes", "stop_TAA", "stop_TAG", "stop_TGA", "frac_TAG", "RF1", "RF2", "RF1_motif", "RF2_motif", "rf_hits"]


def one(acc):
    try:
        contigs = dict(read_fasta(fetch(acc)))
    except Exception as e:  # noqa: BLE001
        return None
    total = sum(len(s) for s in contigs.values())
    finder = pyrodigal.GeneFinder(meta=total < 100000)
    if total >= 100000:
        finder.train(*[s.encode() for s in contigs.values()], translation_table=11)
    stops = Counter()
    prots = {}
    for c, s in contigs.items():
        for g in finder.find_genes(s.encode()):
            nt = g.sequence()
            if not g.partial_end if g.strand == 1 else not g.partial_begin:
                stops[nt[-3:]] += 1
            prots[f"{c}|{g.begin}|{g.end}|{g.strand}"] = g.translate().rstrip("*")
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.faa"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in prots.items()))
        h = Path(td) / "h.tsv"
        subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(RFDB), "-o", str(h), "--threads", "1", "--quiet",
                        "--evalue", "1e-20", "--max-target-seqs", "1", "--outfmt", "6", "qseqid", "sseqid", "pident",
                        "length", "bitscore"], check=True)
        hits = [l.split("\t") for l in open(h)]
    rf1 = rf2 = m1 = m2 = 0
    desc = []
    for qid, sid, pid, ln, bs in hits:
        kind = "RF1" if "|RF1_" in sid else "RF2"
        seq = prots[qid]
        mot1 = bool(re.search(r"P[AVKETSIL]T", seq))
        mot2 = "SPF" in seq
        if kind == "RF1":
            rf1 += 1; m1 += mot1
        else:
            rf2 += 1; m2 += mot2
        desc.append(f"{kind}:{pid}%/{ln}aa{'/PxT' if mot1 else ''}{'/SPF' if mot2 else ''}")
    n = sum(stops.values())
    return dict(acc=acc, bp=total, contigs=len(contigs), genes=len(prots), stop_TAA=stops["TAA"], stop_TAG=stops["TAG"],
                stop_TGA=stops["TGA"], frac_TAG=round(stops["TAG"] / n, 4) if n else "", RF1=rf1, RF2=rf2,
                RF1_motif=m1, RF2_motif=m2, rf_hits=";".join(desc))


if __name__ == "__main__":
    accs = [l.split()[0].replace("GB_", "").replace("RS_", "") for l in open(sys.argv[1]) if l.strip()]
    done = set()
    if OUT.exists():
        done = {l.split("\t")[0] for l in open(OUT)}
    else:
        OUT.write_text("\t".join(COLS) + "\n")
    accs = [a for a in accs if a not in done]
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as p, open(OUT, "a") as f:
        for r in p.imap_unordered(one, accs):
            if r:
                f.write("\t".join(str(r[c]) for c in COLS) + "\n")
                f.flush()
