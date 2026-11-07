"""Translate a genomic region (TAG -> X) and align to a control region's translation (local, BLOSUM62); print alignment
and the control residue/codon opposite each X. Usage: python -I region_translate.py <dir> acc:contig:start:end:strand acc:contig:start:end:strand"""
import sys
from Bio import Align
from Bio.Align import substitution_matrices
R = sys.argv[1]
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
COMP = str.maketrans('ACGT', 'TGCA')
def contig(acc, c):
    out = []; on = False
    for l in open(f'{R}/data/{acc}.fna'):
        if l.startswith('>'):
            if on: break
            on = l[1:].split()[0] == c; continue
        if on: out.append(l.strip())
    return ''.join(out).upper()
def get(spec, tagx):
    acc, c, s, e, st = spec.split(':'); s, e = int(s), int(e); seq = contig(acc, c)[s-1:e]
    if st == '-1': seq = seq.translate(COMP)[::-1]
    cods = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    return ''.join('X' if (tagx and cd == 'TAG') else CODE.get(cd, 'X') for cd in cods), cods
t, tc = get(sys.argv[2], True); o, oc = get(sys.argv[3], False)
al = Align.PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62'); al.open_gap_score = -11; al.extend_gap_score = -1; al.mode = 'local'
a = al.align(t, o)[0]; print(a)
for (ts, te), (os_, oe) in zip(*a.aligned):
    for k in range(te - ts):
        if t[ts + k] == 'X' and tc[ts + k] == 'TAG': print(f'TAG at test codon {ts+k+1}: control {o[os_+k]} ({oc[os_+k]}); test context {t[max(0,ts+k-5):ts+k+6]} control {o[max(0,os_+k-5):os_+k+6]}')
