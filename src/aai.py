"""Pairwise average amino-acid identity (reciprocal best DIAMOND hits, >= 50% query coverage) for a set of genomes.
Genomes listed as acc or acc:table (use 15 for genomes that read TAG as sense, so their genes aren't split).
usage: python src/aai.py <acc[:table]> ...  -> matrix on stdout"""
import subprocess
import sys
import tempfile
from itertools import combinations
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import fetch, read_fasta


def proteins(acc, table, path):
    contigs = dict(read_fasta(fetch(acc)))
    finder = pyrodigal.GeneFinder()
    finder.train(*[s.encode() for s in contigs.values()], translation_table=table)
    with open(path, "w") as f:
        for c, s in contigs.items():
            for g in finder.find_genes(s.encode()):
                f.write(f">{c}|{g.begin}\n{g.translate().rstrip('*')}\n")


def best(q, d, td):
    o = Path(td) / "o.tsv"
    subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(d), "-o", str(o), "--threads", "3", "--quiet",
                    "--evalue", "1e-5", "--max-target-seqs", "1", "--max-hsps", "1", "--query-cover", "50",
                    "--outfmt", "6", "qseqid", "sseqid", "pident"], check=True)
    return {l.split("\t")[0]: (l.split("\t")[1], float(l.split("\t")[2])) for l in open(o)}


def main(specs):
    gs = [(s.split(":")[0], int(s.split(":")[1]) if ":" in s else 11) for s in specs]
    with tempfile.TemporaryDirectory() as td:
        for acc, t in gs:
            p = Path(td) / f"{acc}.faa"
            proteins(acc, t, p)
            subprocess.run(["diamond", "makedb", "--in", str(p), "-d", str(p), "--quiet"], check=True)
        A = {}
        for (a, _), (b, _) in combinations(gs, 2):
            ab = best(Path(td) / f"{a}.faa", Path(td) / f"{b}.faa.dmnd", td)
            ba = best(Path(td) / f"{b}.faa", Path(td) / f"{a}.faa.dmnd", td)
            ids = [pid for q, (s, pid) in ab.items() if ba.get(s, (None,))[0] == q]
            A[a, b] = A[b, a] = (sum(ids) / len(ids) if ids else 0, len(ids))
    names = [a for a, _ in gs]
    print("AAI (%), reciprocal best hits")
    print("\t" + "\t".join(n[4:13] for n in names))
    for a in names:
        print(a[4:13] + "\t" + "\t".join("-" if a == b else f"{A[a, b][0]:.1f}" for b in names))


if __name__ == "__main__":
    main(sys.argv[1:])
