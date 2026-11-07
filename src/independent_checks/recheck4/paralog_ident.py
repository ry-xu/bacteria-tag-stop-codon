"""For given conserved query genes, project every locus (HSP e<1e-10) in each genome onto query coordinates and print the
pairwise % identity between loci over shared query positions (used to see which copy of a duplicated marker in a test
genome is the 'native' one). Usage: python3 -I paralog_ident.py <dir> query_id acc..."""
import sys, itertools
R, Q = sys.argv[1], sys.argv[2]; accs = sys.argv[3:]
loci = {}
for a in accs:
    for l in open(f'{R}/out/{a}.tblastn.tsv'):
        f = l.rstrip('\n').split('\t')
        if f[0] != Q or float(f[8]) >= 1e-10: continue
        s = f[1].split('|')[1] if '|' in f[1] else f[1]
        qs, qseq, sseq = int(f[4]), f[13], f[14]; m = {}; p = qs
        for x, y in zip(qseq, sseq):
            if x != '-':
                if y != '-': m[p] = y
                p += 1
        loci[f'{a[:15]}|{s}:{f[6]}'] = m
ks = list(loci)
for x, y in itertools.combinations(ks, 2):
    sh = set(loci[x]) & set(loci[y])
    if len(sh) < 40: continue
    print(f'{Q}\t{x}\t{y}\tshared={len(sh)}\tident={sum(loci[x][p]==loci[y][p] for p in sh)/len(sh):.3f}')
