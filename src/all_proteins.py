"""Write all pyrodigal proteins of a genome (table 11 and, optionally, a table where one stop is read as X).
usage: python src/all_proteins.py <acc> <out.faa> [TAG|TGA|TAA]"""
import sys
from pathlib import Path
import pyrodigal
sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import fetch, read_fasta, translate, revcomp
acc, out = sys.argv[1], sys.argv[2]
contigs = dict(read_fasta(fetch(acc)))
total = sum(len(s) for s in contigs.values())
finder = pyrodigal.GeneFinder(meta=total < 100000)
if total >= 100000:
    finder.train(*[s.encode() for s in contigs.values()], translation_table=11)
with open(out, "w") as f:
    for c, s in contigs.items():
        for i, g in enumerate(finder.find_genes(s.encode())):
            f.write(f">{c}|{g.begin}|{g.end}|{g.strand}\n{g.translate()}\n")
