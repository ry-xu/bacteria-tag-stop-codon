"""Ribosomal-protein copy census (contamination estimate): for each rpl/rps/rpm query (Bs and Ec), merge HSPs (e<1e-10,
qcov>=0.5) into distinct genomic loci per genome; for genomes with >1 locus print each locus with its contig, TAG count in
the contig (census, out/<acc>.contig_aligned.tsv), and its identity to the best locus of a reference genome (projected
on query coordinates). Usage: python3 -I rp_copies.py <dir> ref_acc acc..."""
import sys, re, collections
R, REF = sys.argv[1], sys.argv[2]; accs = sys.argv[3:]
pat = re.compile(r'^(rp[lsm][A-Z]+)$')
def loci_of(a):
    L = collections.defaultdict(list)
    for l in open(f'{R}/out/{a}.tblastn.tsv'):
        f = l.rstrip('\n').split('\t'); q = f[0]
        if not q.startswith('Bs|'): continue
        gn = q.split('|')[1]
        if not pat.match(gn): continue
        qs, qe, ss, se, ev, qlen = int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[8]), int(f[10])
        if ev >= 1e-10 or (qe - qs + 1) / qlen < 0.5: continue
        s = f[1].split('|')[1] if '|' in f[1] else f[1]
        m = {}; p = qs
        for x, y in zip(f[13], f[14]):
            if x != '-':
                if y != '-': m[p] = y
                p += 1
        L[gn].append((s, min(ss, se), max(ss, se), float(f[9]), m))
    out = {}
    for gn, hs in L.items():
        hs.sort(); mg = []
        for h in hs:
            if mg and mg[-1][0] == h[0] and h[1] <= mg[-1][2]:
                if h[3] > mg[-1][3]: mg[-1] = list(h)
            else: mg.append(list(h))
        out[gn] = mg
    return out
ref = loci_of(REF)
def ident(m1, m2):
    sh = set(m1) & set(m2); return (sum(m1[p] == m2[p] for p in sh) / len(sh)) if len(sh) >= 30 else float('nan')
for a in accs:
    L = loci_of(a)
    tag = {}
    for l in list(open(f'{R}/out/{a}.contig_aligned.tsv'))[1:]:
        f = l.split('\t'); tag[f[0]] = int(f[4])
    multi = {g: v for g, v in L.items() if len(v) > 1}
    print(f'## {a}: Bs ribosomal-protein genes found {len(L)}, with >1 locus: {len(multi)}')
    for g, v in sorted(multi.items()):
        for s, b, e, bits, m in v:
            r = max((ident(m, x[4]) for x in ref.get(g, [])), default=float('nan'))
            print(f'   {g:6s} {s}:{b}-{e} bits={bits:.0f} id_to_{REF[:15]}={r:.2f} contig_TAG={tag.get(s, 0)}')
