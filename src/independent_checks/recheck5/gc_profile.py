"""Per-genome GC and per-contig GC distribution (length-weighted deciles); for T also lists GC of named contigs.
Usage: python3 -I gc_profile.py <dir> acc... [-- contig ...]"""
import sys
R = sys.argv[1]; args = sys.argv[2:]
cont = args[args.index('--') + 1:] if '--' in args else []; accs = args[:args.index('--')] if '--' in args else args
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
for a in accs:
    g = read_fasta(f'{R}/data/{a}.fna'); L = sum(map(len, g.values()))
    gc = {c: (s.count('G') + s.count('C')) / max(1, len(s) - s.count('N')) for c, s in g.items()}
    tot = sum(gc[c] * len(g[c]) for c in g) / L
    xs = sorted(g, key=gc.get); cum = 0; q = []
    for c in xs:
        cum += len(g[c])
        while len(q) < 9 and cum >= L * (len(q) + 1) / 10: q.append(gc[c])
    bins = {}
    for c in g: b = int(gc[c] * 20) / 20; bins[b] = bins.get(b, 0) + len(g[c])
    print(f'{a}: {len(g)} contigs {L/1e6:.2f} Mb GC={tot:.3f} | length-weighted GC deciles ' + ' '.join(f'{x:.2f}' for x in q)
          + ' | Mb per GC bin: ' + ' '.join(f'{b:.2f}:{v/1e6:.2f}' for b, v in sorted(bins.items())))
    for c in cont:
        if c in g: print(f'   {c} len={len(g[c])} GC={gc[c]:.3f}')
