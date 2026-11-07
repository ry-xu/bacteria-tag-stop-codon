"""pyrodigal gene calls (single mode, trained per genome, table 11) for control genomes -> out/genes/<acc>.genes.tsv
(contig, begin, end, strand; 1-based inclusive). Used to define the 'coding' set in tag_meaning.py (control codon must lie
in-frame inside a called control gene). Usage: uv run --no-project --with pyrodigal python -I ctrl_genes.py <dir> acc..."""
import sys, os
import pyrodigal
R = sys.argv[1]; os.makedirs(f'{R}/out/genes', exist_ok=True)
def read_fasta(p):
    d = {}; k = None; buf = []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if k: d[k] = ''.join(buf).upper()
            k = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if k: d[k] = ''.join(buf).upper()
    return d
for acc in sys.argv[2:]:
    g = read_fasta(f'{R}/data/{acc}.fna'); gf = pyrodigal.GeneFinder(); gf.train(*[x.encode() for x in g.values() if len(x) > 20000][:60])
    n = 0
    with open(f'{R}/out/genes/{acc}.genes.tsv', 'w') as o:
        for c, s in g.items():
            for x in gf.find_genes(s.encode()): o.write(f'{c}\t{x.begin}\t{x.end}\t{x.strand}\n'); n += 1
    print(acc, n, 'genes')
