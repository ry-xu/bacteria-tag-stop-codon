"""For marker families with more loci in T than in every CAKUFA01 relative (census tblastn: e<1e-10, qcov>=0.5; Bs and Ec
queries), project each T locus and each locus of the comparison genomes onto query coordinates and report, per T locus, its
best identity to any locus of each comparison genome (shared query positions >= 30), the contig's census TAG count, and
contig length/GC. Usage: python3 -I dup_copies.py <dir> T gene1,gene2,... cmp1,cmp2,..."""
import sys, collections
R, T, genes, cmps = sys.argv[1], sys.argv[2], sys.argv[3].split(','), sys.argv[4].split(',')
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
def loci_of(a, gene):
    L = []
    for l in open(f'{R}/out/{a}.tblastn.tsv'):
        f = l.rstrip('\n').split('\t'); q = f[0]
        if not (q.startswith('Ec|') or q.startswith('Bs|')) or q.split('|')[1] != gene: continue
        qs, qe, ss, se, ev, qlen = int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[8]), int(f[10])
        if ev >= 1e-10 or (qe - qs + 1) / qlen < 0.5: continue
        s = f[1].split('|')[1] if '|' in f[1] else f[1]
        m = {}; p = qs
        for x, y in zip(f[13], f[14]):
            if x != '-':
                if y != '-': m[(q[:2], p)] = y
                p += 1
        L.append((s, min(ss, se), max(ss, se), float(f[9]), q[:2], m))
    L.sort(); mg = []
    for h in L:
        if mg and mg[-1][0] == h[0] and h[1] <= mg[-1][2]:
            mg[-1][2] = max(mg[-1][2], h[2]); mg[-1][5].update(h[5])
        else: mg.append([h[0], h[1], h[2], h[3], h[4], dict(h[5])])
    return mg
def ident(m1, m2):
    sh = set(m1) & set(m2); return (sum(m1[p] == m2[p] for p in sh) / len(sh), len(sh)) if len(sh) >= 30 else (float('nan'), len(sh))
g = read_fasta(f'{R}/data/{T}.fna')
tag = {}
for l in list(open(f'{R}/out/{T}.contig_aligned.tsv'))[1:]:
    f = l.split('\t'); tag[f[0]] = int(f[4])
for gene in genes:
    tl = loci_of(T, gene); cl = {c: loci_of(c, gene) for c in cmps}
    print(f'## {gene}: T loci {len(tl)}; ' + ' '.join(f'{c}:{len(v)}' for c, v in cl.items()))
    for s, a, b, bits, src, m in tl:
        cells = []
        for c in cmps:
            best = max((ident(m, x[5]) for x in cl[c]), key=lambda t: (t[0] if t[0] == t[0] else -1), default=(float('nan'), 0))
            cells.append(f'{c[4:13]}={best[0]:.2f}')
        seq = g[s]; gc = (seq.count('G') + seq.count('C')) / len(seq)
        print(f'   {s}:{a}-{b} contig_len={len(seq)} GC={gc:.3f} contig_census_TAG={tag.get(s, 0)} | ' + ' '.join(cells))
