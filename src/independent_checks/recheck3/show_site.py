"""Show +-N codons around a column from an out/cols table: test translation (TAG->X) and the control's translation in its
aligned frame, plus downstream distance to the next stop in both. Usage: python -I show_site.py <dir> T O N [filter: star|all]"""
import sys
R, T, O, N = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]); filt = sys.argv[5] if len(sys.argv) > 5 else 'star'
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
COMP = str.maketrans('ACGT', 'TGCA')
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
gT = read_fasta(f'{R}/data/{T}.fna'); gO = read_fasta(f'{R}/data/{O}.fna')
def window(g, c, g0, strand, n):
    s = g[c]
    if strand == 1: cods = [s[g0 + 3*k: g0 + 3*k + 3] for k in range(-n, n + 1) if 0 <= g0 + 3*k and g0 + 3*k + 3 <= len(s)]
    else: cods = [s[g0 - 3*k: g0 - 3*k + 3].translate(COMP)[::-1] for k in range(-n, n + 1) if 0 <= g0 - 3*k and g0 - 3*k + 3 <= len(s)]
    return ''.join(('x' if (cd == 'TAG') else CODE.get(cd, '?')) if True else '' for cd in cods)
for l in list(open(f'{R}/out/cols/{T}__{O}.tsv'))[1:]:
    r = l.rstrip('\n').split('\t')
    if r[3] != 'TAG' or not (float(r[11]) >= 0.5 and int(r[12]) >= 10 and int(r[13]) == 0): continue
    if filt == 'star' and r[7] != '*': continue
    c, g0, st = r[0], int(r[1]), int(r[2])
    # O: codon start O_pos is 1-based start on + strand (if + frame) or end (if - frame); determine orientation by codon match
    oc, op, ocod = r[4], int(r[5]), r[6]
    s = gO[oc]
    if s[op-1:op+2] == ocod: ow = window(gO, oc, op - 1, 1, N)
    else: ow = window(gO, oc, op - 3, -1, N)
    print(f'{c}:{g0}:{st} O={oc}:{op} {ocod} pid={r[9]}\n   T {window(gT, c, g0, st, N)}\n   O {ow}')
