"""Chimera check without coverage: per test-genome contig, GC, length, genome-wide TAG density (TAG codons among all
aligned test codons in out/cols/<T>__<O>.tsv, union over the given other genomes) and per-contig amino-acid identity to each
other genome (length-weighted pident of the best HSP per test ORF, from out/cols/<T>__<O>.tblastn.tsv).
Groups contigs by identity-to-sister quantile and by GC tercile and reports TAG density per group; lists contigs with
>=300 aligned codons and zero TAG (expected count given the genome-wide rate), and a Poisson dispersion test.
Usage: python3 -I contig_origin.py <dir> T sister O1,O2,..."""
import sys, collections, math
R, T, SIS, others = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(',')
allo = [SIS] + others
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
g = read_fasta(f'{R}/data/{T}.fna')
# ORF -> contig map from the ORF fasta written by codon_columns.py is not kept per contig; rebuild from cols tsv instead:
cod = {}   # (contig,g0,strand) -> codon ; union over others
for O in allo:
    for l in list(open(f'{R}/out/cols/{T}__{O}.tsv'))[1:]:
        f = l.split('\t')
        if float(f[11]) >= 0.4 and int(f[12]) >= 5: cod[(f[0], f[1], f[2])] = f[3]
n_al = collections.Counter(k[0] for k in cod); n_tag = collections.Counter(k[0] for k, c in cod.items() if c == 'TAG')
# per-contig identity: aligned columns per contig per O with identity
idn = {O: collections.defaultdict(lambda: [0, 0]) for O in allo}
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
for O in allo:
    for l in list(open(f'{R}/out/cols/{T}__{O}.tsv'))[1:]:
        f = l.split('\t')
        if int(f[12]) < 5: continue
        a = 'X' if f[3] == 'TAG' else CODE.get(f[3], 'X')
        idn[O][f[0]][0] += 1; idn[O][f[0]][1] += (a == f[7])
rate = sum(n_tag.values()) / sum(n_al.values())
gc = {c: (s.count('G') + s.count('C')) / len(s) for c, s in g.items()}
rows = []
for c in g:
    if n_al[c] < 100: continue
    ids = {O: (idn[O][c][1] / idn[O][c][0] if idn[O][c][0] >= 50 else float('nan')) for O in allo}
    rows.append((c, len(g[c]), gc[c], n_al[c], n_tag[c], ids))
print(f'## {T}: contigs with >=100 aligned codons: {len(rows)} of {len(g)}; genome-wide TAG density {rate*1e5:.0f}/100k aligned codons '
      f'({sum(n_tag.values())} TAG / {sum(n_al.values())})')
chi = sum((r[4] - r[3] * rate) ** 2 / (r[3] * rate) for r in rows)
print(f'   dispersion chi2={chi:.1f} df={len(rows)-1}')
def grp(lab, key):
    rs = sorted(rows, key=key); k = len(rs) // 3
    for name, part in (('low', rs[:k]), ('mid', rs[k:2*k]), ('high', rs[2*k:])):
        A = sum(r[3] for r in part); Tg = sum(r[4] for r in part)
        print(f'   {lab} tercile {name:4s}: contigs={len(part)} range={key(part[0]):.3f}-{key(part[-1]):.3f} aligned={A} TAG={Tg} density={Tg/A*1e5:.0f}/100k')
grp('GC', lambda r: r[2])
for O in allo:
    grp(f'aa-identity to {O}', lambda r, O=O: (r[5][O] if r[5][O] == r[5][O] else -1))
z = [r for r in rows if r[4] == 0 and r[3] * rate >= 3]
print(f'   contigs with zero TAG where >=3 expected: {len(z)}; P(0|exp) product = {math.exp(-sum(r[3]*rate for r in z)):.2g}')
for r in sorted(z, key=lambda r: -r[3]):
    print(f'     {r[0]} len={r[1]} GC={r[2]:.3f} aligned={r[3]} exp_TAG={r[3]*rate:.1f} ' + ' '.join(f'id_{O[:13]}={r[5][O]:.2f}' for O in allo))
print('   top TAG-dense contigs:')
for r in sorted(rows, key=lambda r: -(r[4] - r[3] * rate) / math.sqrt(r[3] * rate))[:6]:
    print(f'     {r[0]} len={r[1]} GC={r[2]:.3f} aligned={r[3]} TAG={r[4]} exp={r[3]*rate:.1f} ' + ' '.join(f'id_{O[:13]}={r[5][O]:.2f}' for O in allo))
with open(f'{R}/out/{T}.contig_origin.tsv', 'w') as o:
    o.write('contig\tlen\tGC\taligned_codons\tTAG\t' + '\t'.join(f'aaid_{O}' for O in allo) + '\n')
    for r in rows: o.write(f'{r[0]}\t{r[1]}\t{r[2]:.3f}\t{r[3]}\t{r[4]}\t' + '\t'.join(f'{r[5][O]:.3f}' for O in allo) + '\n')
