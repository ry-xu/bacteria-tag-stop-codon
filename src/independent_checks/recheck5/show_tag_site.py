"""Show the T tblastn HSP (from out/cols/<T>__<O>.tblastn.tsv) that covers a given T codon, with T's translated ORF
(TAG = X) aligned to O's translation (std code, '*' = stop), +-W columns around the site, and the T ORF coordinates.
Usage: python3 -I show_tag_site.py <dir> T O contig g0 strand [W]"""
import sys
R, T, O, c, g0, st = sys.argv[1:7]; g0 = int(g0); st = int(st); W = int(sys.argv[7]) if len(sys.argv) > 7 else 30
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
gT = read_fasta(f'{R}/data/{T}.fna'); L = len(gT[c])
# rebuild ORF index exactly as codon_columns.py
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
COMP = str.maketrans('ACGT', 'TGCA'); rc = lambda s: s.translate(COMP)[::-1]
orfs = {}; n = 0
for cc, s in gT.items():
    for strand, seq in ((1, s), (-1, rc(s))):
        for fr in range(3):
            plen = 0; start = fr
            for p in range(fr, len(s) - 2, 3):
                cod = seq[p:p+3]; aa = 'X' if cod == 'TAG' else CODE.get(cod, 'X')
                if aa == '*':
                    if plen >= 100: orfs[f'o{n}'] = (cc, strand, start, plen); n += 1
                    plen = 0; start = p + 3
                else: plen += 1
            if plen >= 100: orfs[f'o{n}'] = (cc, strand, start, plen); n += 1
o_ = st == 1 and g0 or (L - g0 - 3)
hit = [k for k, (cc, s_, a, pl) in orfs.items() if cc == c and s_ == st and a <= o_ < a + 3 * pl and (o_ - a) % 3 == 0]
name = hit[0]; cc, s_, a, pl = orfs[name]; qpos = (o_ - a) // 3 + 1
print(f'site {c}:{g0}:{st} -> ORF {name} ({pl} aa, starts at oriented offset {a}); query position {qpos}')
best = None
for l in open(f'{R}/out/cols/{T}__{O}.tblastn.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] != name or not (int(f[4]) <= qpos <= int(f[5])): continue
    if best is None or float(f[9]) > float(best[9]): best = f
f = best; qs = int(f[4]); qseq, sseq = f[11], f[12]
p = qs; k0 = None
for k, x in enumerate(qseq):
    if x != '-':
        if p == qpos: k0 = k; break
        p += 1
lo, hi = max(0, k0 - W), min(len(qseq), k0 + W + 1)
print(f'HSP vs {O} {f[1]} q{f[4]}-{f[5]} s{f[6]}-{f[7]} pid={f[2]} bits={f[9]}')
print('T  ' + qseq[lo:k0] + '[' + qseq[k0] + ']' + qseq[k0+1:hi])
print('O  ' + sseq[lo:k0] + '[' + sseq[k0] + ']' + sseq[k0+1:hi])
