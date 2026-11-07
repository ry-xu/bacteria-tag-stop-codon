"""Per-contig spread of in-alignment TAG codons: Poisson expectation from aligned codons per contig,
dispersion index, GC of TAG-bearing vs other contigs. Usage: python -I contig_spread.py <recheck_dir> acc..."""
import sys, math, collections
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
for acc in sys.argv[2:]:
    g = read_fasta(f'{R}/data/{acc}.fna')
    rows = [l.rstrip('\n').split('\t') for l in open(f'{R}/out/{acc}.contig_aligned.tsv')][1:]
    N = sum(int(r[2]) for r in rows); T = sum(int(r[4]) for r in rows); rate = T / N
    chi = sum((int(r[4]) - int(r[2]) * rate) ** 2 / (int(r[2]) * rate) for r in rows)
    withtag = [r for r in rows if int(r[4]) > 0]
    gc = lambda c: (g[c].count('G') + g[c].count('C')) / len(g[c])
    def wgc(cs):
        L = sum(len(g[c]) for c in cs); return sum(gc(c) * len(g[c]) for c in cs) / L
    allgc = wgc(list(g)); tg = wgc([r[0] for r in withtag]); og = wgc([c for c in g if c not in {r[0] for r in withtag}])
    top = max(rows, key=lambda r: int(r[4]))
    # Poisson tail for the top contig
    lam = int(top[2]) * rate; k = int(top[4])
    p = 1 - sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k))
    print(f'{acc}: contigs_with_alignments={len(rows)} contigs_with_TAG={len(withtag)} TAG={T} aligned_codons={N} '
          f'rate={rate*1e5:.1f}/100k chi2={chi:.1f} df={len(rows)-1} | max contig {top[0]} TAG={k} exp={lam:.2f} P(>=k)={p:.2g} (x{len(rows)} contigs) '
          f'| GC all={allgc:.3f} TAG-contigs={tg:.3f} other={og:.3f} | TAG-contigs carry {sum(int(r[2]) for r in withtag)/N:.0%} of aligned codons')
