"""(Adapted from recheck/scripts; loci = this run's HMM prfA calls.) Extract stop-to-stop ORFs (>=100 aa, std code, trimmed to first M) within +-FLANK nt of each
control genome's prfA locus; write them as proteins named by order along the control contig.
Usage: python -I prfA_neighborhood.py <recheck_dir>"""
import sys
R = sys.argv[1]; FLANK = 15000
loci = {'GCA_934796175.1': ('CAKVOV010000023.1', 30989, 32167, 1),
        'GCA_004556155.1': ('SFIM01000024.1', 3755, 4933, -1),
        'GCA_963624325.1': ('CAVFPH010000002.1', 73650, 74726, 1)}
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
out = open(f'{R}/data/prfA_neighbors.faa', 'w')
for acc, (c, s, e, st) in loci.items():
    seq = read_fasta(f'{R}/data/{acc}.fna')[c]
    lo, hi = max(0, s - 1 - FLANK), min(len(seq), e + FLANK)
    reg = seq[lo:hi]; L = len(reg); orfs = []
    for strand, sq in ((1, reg), (-1, reg.translate(COMP)[::-1])):
        for fr in range(3):
            prot = ''; start = fr
            for p in range(fr, L - 2, 3):
                aa = CODE.get(sq[p:p+3], 'X')
                if aa == '*':
                    m = prot.find('M')
                    if m >= 0 and len(prot) - m >= 100:
                        a0 = start + 3 * m; a1 = p + 3
                        g0, g1 = (lo + a0 + 1, lo + a1) if strand == 1 else (lo + L - a1 + 1, lo + L - a0)
                        orfs.append((g0, g1, strand, prot[m:]))
                    prot = ''; start = p + 3
                else: prot += aa
    orfs.sort()
    for g0, g1, strand, prot in orfs:
        tag = "prfA" if (g0 <= e and g1 >= s and min(g1, e) - max(g0, s) > 0.5 * (e - s)) else ""
        out.write(f'>{acc}|{c}|{g0}-{g1}|{strand}|{tag}\n{prot}\n')
out.close()
