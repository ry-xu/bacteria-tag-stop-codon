"""Global identity of the CTA tRNA vs every tRNA called in the five MAGs (top hits), plus flanks of the
CTA tRNA locus in both test genomes. Usage: uv run --with biopython python -I trna_compare.py <dir>"""
import sys
from Bio import Align
R = sys.argv[1]
def rf(p):
    d = {}; k = None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = ''
        else: d[k] += l.upper()
    return d
t = rf(f'{R}/out/trna/all_trna.fa'); sup = rf(f'{R}/out/trna/sup_cta.fa')
s = list(sup.values())[0]
al = Align.PairwiseAligner(); al.mode = 'global'; al.match_score = 1; al.mismatch_score = -1; al.open_gap_score = -2; al.extend_gap_score = -1
res = []
for k, v in t.items():
    if 'CTA' in k: continue
    a = al.align(s, v)[0]; m = sum(1 for (x, y) in zip(*a) if x == y and x != '-')
    res.append((m / max(len(s), len(v)), k, len(v)))
res.sort(reverse=True)
for r in res[:12]: print(f'{r[0]:.3f}\t{r[1]}\t{r[2]}')
# flanks
def rg(p):
    d = {}; k = None; buf = []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if k: d[k] = ''.join(buf).upper()
            k = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if k: d[k] = ''.join(buf).upper()
    return d
C = str.maketrans('ACGT', 'TGCA')
loc = {'GCA_024404275.1': ('JAKSQG010000048.1', 5836, 5911, -1), 'GCA_045599485.1': ('JAZRWA010000050.1', 8748, 8823, 1)}
fl = {}
for acc, (c, a, b, st) in loc.items():
    g = rg(f'{R}/data/{acc}.fna')[c]; F = 300
    seg = g[max(0, a - 1 - F): b + F]
    if st == -1: seg = seg.translate(C)[::-1]
    fl[acc] = seg; print(acc, c, 'contig_len', len(g), 'locus', a, b, st)
x, y = fl.values()
up_x, up_y = x[:300], y[:300]; dn_x, dn_y = x[-300:], y[-300:]
for nm, p, q in (('upstream300', up_x, up_y), ('downstream300', dn_x, dn_y)):
    a = al.align(p, q)[0]; m = sum(1 for (u, v) in zip(*a) if u == v and u != '-'); print(nm, 'identity', round(m / 300, 3))
open(f'{R}/out/trna/sup_cta_flanks.fa', 'w').write(''.join(f'>{k}\n{v}\n' for k, v in fl.items()))
