"""Closest genomes to a query genome among a list: proteins of every listed genome pooled into one DIAMOND database,
query proteins searched against it; per target genome, AAI-like score = mean identity of the query's best hit in that
genome (>= 50% query coverage), and the number of such hits.
usage: python src/nearest_genomes.py <query_acc[:table]> <list: acc<TAB>label> [top=15]"""
import subprocess
import sys
import tempfile
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import fetch, read_fasta


def proteins(task):
    acc, table = task
    try:
        contigs = [s.encode() for _, s in read_fasta(fetch(acc))]
    except Exception as e:  # noqa: BLE001
        return acc, f"ERROR {e}"
    f = pyrodigal.GeneFinder()
    f.train(*contigs, translation_table=table)
    out = []
    for i, c in enumerate(contigs):
        for j, g in enumerate(f.find_genes(c)):
            out.append((f"{acc}#{i}_{j}", g.translate(translation_table=table).rstrip("*")))
    return acc, out


def main():
    q = sys.argv[1].split(":")
    qacc, qtab = q[0], int(q[1]) if len(q) > 1 else 11
    label = dict(l.rstrip("\n").split("\t")[:2] for l in open(sys.argv[2]) if l.strip())
    top = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    with Pool(3) as p:
        res = dict(p.map(proteins, [(qacc, qtab)] + [(a, 11) for a in label if a != qacc]))
    with tempfile.TemporaryDirectory() as td:
        db, qf, o = Path(td) / "db.faa", Path(td) / "q.faa", Path(td) / "o.tsv"
        with open(db, "w") as f:
            for a, ps in res.items():
                if a != qacc and not isinstance(ps, str):
                    f.write("".join(f">{i}\n{s}\n" for i, s in ps))
        qf.write_text("".join(f">{i}\n{s}\n" for i, s in res[qacc]))
        subprocess.run(["diamond", "makedb", "--in", db, "-d", Path(td) / "db", "--quiet"], check=True)
        subprocess.run(["diamond", "blastp", "-q", qf, "-d", Path(td) / "db", "-o", o, "--quiet", "--threads", "4",
                        "--max-target-seqs", "2000", "--evalue", "1e-10", "--query-cover", "50",
                        "--outfmt", "6", "qseqid", "sseqid", "pident"], check=True)
        best = defaultdict(dict)
        for line in open(o):
            qi, si, pid = line.split("\t")
            g = si.split("#")[0]
            if qi not in best[g] or float(pid) > best[g][qi]:
                best[g][qi] = float(pid)
    rows = sorted(((sum(v.values()) / len(v), len(v), g) for g, v in best.items() if len(v) >= 100), reverse=True)
    print(f"{qacc}: {len(res[qacc])} proteins; {len(rows)} genomes with >= 100 shared proteins")
    for aai, n, g in rows[:top]:
        print(f"  {g}  AAI~{aai:.1f}  n={n}  {label.get(g, '')}")


if __name__ == "__main__":
    main()
