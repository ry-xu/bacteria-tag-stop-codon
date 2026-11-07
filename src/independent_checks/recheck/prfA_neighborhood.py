"""Extract stop-to-stop ORFs (>=100 aa, std code, trimmed to first M) within +-FLANK nt of each
control genome's prfA locus; write them as proteins named by order along the control contig.
Usage: python -I prfA_neighborhood.py <recheck_dir>"""
import sys
R = sys.argv[1]; FLANK = 15000
loci = {'GCA_015061915.1': ('SUVR01000013.1', 51355, 52527, 1),
        'GCA_934644965.1': ('CAKUDG010000004.1', 36255, 37472, 1),
        'GCA_045589595.1': ('JAZNPK010000002.1', 252748, 253905, 1)}
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
        tag = 'prfA' if abs(g0 - s) < 200 or abs(g1 - e) < 200 else ''
        out.write(f'>{acc}|{c}|{g0}-{g1}|{strand}|{tag}\n{prot}\n')
out.close()
