"""Gene order in genomic windows: pyrodigal (trained on the whole genome, table 11; for TAG-test genomes also a run
with TAG not a stop is impossible in pyrodigal, so genes are std-code calls) + DIAMOND top Swiss-Prot hit.
Writes the window DNA to out/rf/win_<label>.fna and genes to out/rf/win_<label>.faa.
Usage: uv run --no-project --with pyrodigal python -I locus_genes.py <dir> label:acc:contig:start:end ..."""
import sys, subprocess
import pyrodigal
R = sys.argv[1]
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
cache = {}
for spec in sys.argv[2:]:
    lab, acc, c, s, e = spec.split(':'); s, e = int(s), int(e)
    if acc not in cache:
        g = read_fasta(f'{R}/data/{acc}.fna'); gf = pyrodigal.GeneFinder(); gf.train(*[x.encode() for x in g.values() if len(x) > 20000][:50]); cache[acc] = (g, gf)
    g, gf = cache[acc]; seq = g[c]; e = min(e, len(seq))
    genes = [x for x in gf.find_genes(seq.encode()) if x.end >= s and x.begin <= e]
    with open(f'{R}/out/rf/win_{lab}.fna', 'w') as o: o.write(f'>{lab}|{acc}|{c}|{s}-{e}\n{seq[s-1:e]}\n')
    fa = f'{R}/out/rf/win_{lab}.faa'
    with open(fa, 'w') as o:
        for x in genes: o.write(f'>{x.begin}-{x.end}|{x.strand:+d}\n{x.translate().rstrip("*")}\n')
    out = subprocess.run(['diamond', 'blastp', '-q', fa, '-d', f'{R}/data/sprot', '--threads', '2', '--quiet', '-k', '1', '-e', '1e-5',
                          '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'stitle'], capture_output=True, text=True, check=True).stdout
    ann = {}
    for l in out.splitlines():
        q, sid, pid, st = l.split('\t')
        ann.setdefault(q, f'{st.split(" OS=")[0].split(" ", 1)[1][:50]} [{sid.split("|")[2]} {float(pid):.0f}%]')
    print(f'== {lab} {acc} {c} ({len(seq)} nt) window {s}-{e}')
    prev = None
    for x in genes:
        k = f'{x.begin}-{x.end}|{x.strand:+d}'
        gap = '' if prev is None else f'gap {x.begin - prev - 1:+d}'
        print(f'   {x.begin:>7}-{x.end:<7} {x.strand:+d} {len(x.translate())-1:>4}aa stop={seq[x.end-3:x.end] if x.strand==1 else ""}{"".join({"A":"T","C":"G","G":"C","T":"A"}[b] for b in seq[x.begin-1:x.begin+2][::-1]) if x.strand==-1 else ""} {gap:>9}  {ann.get(k, "-")}')
        prev = x.end
