"""Pyrodigal (single mode, table 11) genes for each genome in data/*.fna -> out/genes/<g>.faa, <g>.ffn, <g>.tsv
Protein ids: <genome>|<contig>|<begin>|<end>|<strand>|<stopcodon>"""
import sys
from pathlib import Path
import pyrodigal

D = Path(__file__).resolve().parent.parent
def read_fasta(p):
    name, seq = None, []
    for line in open(p):
        if line[0] == ">":
            if name: yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip().upper())
    if name: yield name, "".join(seq)
(D / "out" / "genes").mkdir(parents=True, exist_ok=True)
for g in sys.argv[1:]:
    contigs = dict(read_fasta(D / "data" / f"{g}.fna"))
    gf = pyrodigal.GeneFinder()
    gf.train(*[s.encode() for s in contigs.values()], translation_table=11)
    with open(D / "out" / "genes" / f"{g}.faa", "w") as fa, open(D / "out" / "genes" / f"{g}.ffn", "w") as fn, \
         open(D / "out" / "genes" / f"{g}.tsv", "w") as tsv:
        tsv.write("id\tcontig\tbegin\tend\tstrand\tstop\tpartial\n")
        n = 0
        for c, s in contigs.items():
            for gene in gf.find_genes(s.encode()):
                nt = gene.sequence()
                stop = nt[-3:] if not gene.partial_end else "NA"
                gid = f"{g}|{c}|{gene.begin}|{gene.end}|{gene.strand}|{stop}"
                prot = gene.translate(translation_table=11).rstrip("*")
                fa.write(f">{gid}\n{prot}\n"); fn.write(f">{gid}\n{nt}\n")
                tsv.write(f"{gid}\t{c}\t{gene.begin}\t{gene.end}\t{gene.strand}\t{stop}\t{int(gene.partial_begin)}{int(gene.partial_end)}\n")
                n += 1
    print(g, len(contigs), n, flush=True)
