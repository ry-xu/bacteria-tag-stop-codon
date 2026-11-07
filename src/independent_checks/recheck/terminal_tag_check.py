"""For HSPs reaching within 10 aa of the query C-terminus, find the first in-frame stop after the HSP end;
when that is TAG, report how far the next TAA/TGA is. Usage: python -I terminal_tag_check.py <recheck_dir> acc..."""
import sys
R = sys.argv[1]; COMP = str.maketrans('ACGT', 'TGCA')
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
    g = read_fasta(f'{R}/data/{acc}.fna'); seen = set()
    for line in open(f'{R}/out/{acc}.tblastn.tsv'):
        f = line.rstrip('\n').split('\t'); q = f[0]
        if not (q.startswith('Ec|') or q.startswith('Bs|')): continue
        s = f[1].split('|')[1] if '|' in f[1] else f[1]
        qs, qe, ss, se, ev, qlen = int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[8]), int(f[10])
        if ev >= 1e-10 or (qe - qs + 1) / qlen < 0.5 or qe < qlen - 10: continue
        seq = g[s]; strand = 1 if se > ss else -1
        p = se + 1 if strand == 1 else se - 1   # first nt after HSP
        cods = []
        for k in range(300):
            if strand == 1:
                if p + 2 > len(seq): break
                c = seq[p-1:p+2]; p += 3
            else:
                if p - 2 < 1: break
                c = seq[p-3:p].translate(COMP)[::-1]; p -= 3
            cods.append(c)
        first = next((i for i, c in enumerate(cods) if c in ('TAA', 'TAG', 'TGA')), None)
        if first is None or cods[first] != 'TAG': continue
        key = (s, strand, first if strand == 1 else -first, se)
        nxt = next((i for i, c in enumerate(cods) if c in ('TAA', 'TGA')), None)
        k2 = (s, se + (3 * first if strand == 1 else -3 * first))
        if k2 in seen: continue
        seen.add(k2)
        print(acc, q, s, f'HSP_end={se}', f'q{qs}-{qe}/{qlen}', f'TAG {first} codons after HSP end;',
              f'next TAA/TGA {"none in 300 codons" if nxt is None else str(nxt - first) + " codons after the TAG"}; TAGs before it: {cods[:nxt if nxt else 300].count("TAG")}')
