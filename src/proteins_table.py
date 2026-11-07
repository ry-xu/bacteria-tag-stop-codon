"""All pyrodigal proteins of a genome under a given NCBI translation table (default 11).
usage: python src/proteins_table.py <acc> <out.faa> [table]"""
import sys
from pathlib import Path
import pyrodigal
sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import fetch, read_fasta
acc, out = sys.argv[1], sys.argv[2]
table = int(sys.argv[3]) if len(sys.argv) > 3 else 11
contigs = dict(read_fasta(fetch(acc)))
total = sum(len(s) for s in contigs.values())
finder = pyrodigal.GeneFinder(meta=total < 100000)
if total >= 100000:
    finder.train(*[s.encode() for s in contigs.values()], translation_table=table)
n = cd = 0
with open(out, "w") as f:
    for c, s in contigs.items():
        for g in finder.find_genes(s.encode()):
            n += 1
            cd += g.end - g.begin + 1
            f.write(f">{c}|{g.begin}|{g.end}|{g.strand}\n{g.translate()}\n")
print(f"{acc} table {table}: {n} genes, coding density {cd / total:.3f}")
