"""What does TAG encode? Uses out/cols/<T>__<O>.tsv (codon_columns.py: test ORFs read with TAG->X, tblastn vs O's raw DNA,
one column per test codon, frame-competition + ungapped +-3 filter).
Sets: 'census' = test codon inside a frame-matched conserved-protein HSP (step-1 set);
      'genome' = any column with win10_ident >= 0.4 and >= 10 columns from the HSP edge;
      'coding' = win10_ident >= 0.5, >= 10 columns from the edge, and the control codon lies in-frame inside a
                 pyrodigal-called control gene (ctrl_genes.py) -> excludes non-coding/antisense alignments;
      'strict' = win10_ident >= 0.5, >= 10 columns from the edge, and no stop in the control's translated HSP within +-20
                 columns (excludes alignments to control non-coding / shadow frames).
Per test genome and its controls:
 A. TAG sites: control residue counts per control; per-site majority over controls.
 B. Calibration: same statistic for test GAG, GAA, CAG, CAA, other codons.
 C. Likelihood: profile P(control residue | test aa A) from all test sense codons of A in the same set (pooled controls,
    add-0.5 smoothing); background profile = control residues at all aligned columns. Log-likelihood of TAG-site
    residues under each A and under background, (i) pooled over controls and (ii) one residue per site (per-site majority,
    ties -> first control) to avoid pseudo-replication; bootstrap over sites for LL(Glu)-LL(Gln) and LL(Glu)-LL(bg).
 D. Mixture: EM for the fraction of TAG sites drawn from each amino-acid profile (+ background), per-site residues.
Usage: python3 -I tag_meaning.py <dir> T ctrl1,ctrl2,..."""
import sys, collections, math, random
R, T, ctrls = sys.argv[1], sys.argv[2], sys.argv[3].split(','); C = f'{R}/out/cols'
random.seed(7)
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
AA = 'ACDEFGHIKLMNPQRSTVWY'
def load(T, O):
    rows = [l.rstrip('\n').split('\t') + [O] for l in open(f'{C}/{T}__{O}.tsv')][1:]
    return {(r[0], r[1], r[2]): r for r in rows}
def sel(r, mode):
    if mode == 'census': return r[10] == '1'
    if mode == 'coding': return float(r[11]) >= 0.5 and int(r[12]) >= 10 and (r[0], r[1], r[2]) in coding_ok[r[-1]]
    if mode == 'genome': return float(r[11]) >= 0.4 and int(r[12]) >= 10
    return float(r[11]) >= 0.5 and int(r[12]) >= 10 and int(r[13]) == 0   # strict: control frame is open (no stop) +-20 columns
def fmt(cnt, k=10):
    n = sum(cnt.values()); return f'n={n} ' + ' '.join(f'{a}:{c}' for a, c in cnt.most_common(k))
tabs = {O: load(T, O) for O in ctrls}
# 'coding' set: control codon lies in-frame inside a pyrodigal gene of the control (out/genes/<O>.genes.tsv)
import bisect
def read_contigs(acc):
    d = {}; k = None; buf = []
    for line in open(f'{R}/data/{acc}.fna'):
        line = line.rstrip()
        if line.startswith('>'):
            if k: d[k] = ''.join(buf).upper()
            k = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if k: d[k] = ''.join(buf).upper()
    return d
COMPL = str.maketrans('ACGT', 'TGCA')
coding_ok = {}
for O in ctrls:
    genes = {}
    for l in open(f'{R}/out/genes/{O}.genes.tsv'):
        c, b_, e_, st = l.split('\t'); genes.setdefault(c, []).append((int(b_), int(e_), int(st)))
    gO = read_contigs(O)
    def in_gene(c, a, z, strand):   # codon at 1-based a..z on given strand
        for b_, e_, st in genes.get(c, ()):
            if st != strand or a < b_ or z > e_: continue
            if (strand == 1 and (a - b_) % 3 == 0) or (strand == -1 and (e_ - z) % 3 == 0): return True
        return False
    ok = set()
    for k, r in tabs[O].items():
        oc, op, ocod = r[4], int(r[5]), r[6]; s_ = gO[oc]
        if s_[op-1:op+2] == ocod and in_gene(oc, op, op + 2, 1): ok.add(k)
        elif op >= 3 and s_[op-3:op].translate(COMPL)[::-1] == ocod and in_gene(oc, op - 2, op, -1): ok.add(k)
    coding_ok[O] = ok
for mode in ('census', 'coding', 'strict', 'genome'):
    print(f'\n===== {T} vs {",".join(ctrls)}  set={mode}')
    sites = sorted(set().union(*[{k for k, r in tabs[O].items() if r[3] == 'TAG' and sel(r, mode)} for O in ctrls]))
    print(f'A. TAG sites with >=1 control aligned: {len(sites)}')
    for O in ctrls:
        cnt = collections.Counter(tabs[O][k][7] for k in sites if k in tabs[O] and sel(tabs[O][k], mode))
        print(f'   vs {O}: {fmt(cnt, 12)}')
    persite = {}
    cons = collections.Counter()
    for k in sites:
        res = [tabs[O][k][7] for O in ctrls if k in tabs[O] and sel(tabs[O][k], mode)]
        top, n = collections.Counter(res).most_common(1)[0]
        cons[top if n > len(res) / 2 else 'mixed'] += 1
        persite[k] = top
    print(f'   per-site majority over controls: {fmt(cons, 14)}')
    print(f'   control codons at TAG sites (pooled): {fmt(collections.Counter(tabs[O][k][6] for O in ctrls for k in sites if k in tabs[O] and sel(tabs[O][k], mode)), 12)}')
    print('B. calibration (pooled over controls; share of control residue):')
    for cods in (('GAG',), ('GAA',), ('CAG',), ('CAA',), ('GAT', 'GAC'), ('AAG', 'AAA'), ('TAT', 'TAC'), ('TGG',), ('TTG',), ('TCG',), ('AGG',), ('TTA',)):
        cnt = collections.Counter(r[7] for O in ctrls for r in tabs[O].values() if r[3] in cods and sel(r, mode))
        n = sum(cnt.values())
        print(f'   test {"/".join(cods):8s}({CODE[cods[0]]}): n={n:6d} ' + ' '.join(f'{a}:{cnt[a]/n:.1%}' for a, _ in cnt.most_common(5)))
    cnt = collections.Counter(r[7] for O in ctrls for r in tabs[O].values() if r[3] == 'TAG' and sel(r, mode))
    n = sum(cnt.values())
    print(f'   test TAG     (?): n={n:6d} ' + ' '.join(f'{a}:{cnt[a]/n:.1%}' for a, _ in cnt.most_common(8)))
    # C. likelihood
    prof = {A: collections.Counter() for A in AA}; bgc = collections.Counter()
    for O in ctrls:
        for r in tabs[O].values():
            if not sel(r, mode) or r[7] not in AA: continue
            a = CODE.get(r[3])
            if a in prof: prof[a][r[7]] += 1
            if r[3] != 'TAG': bgc[r[7]] += 1
    prof['bg'] = bgc
    def lp(A, x):
        n = sum(prof[A].values()); return math.log((prof[A][x] + 0.5) / (n + 10))
    obs_pool = [r[7] for O in ctrls for k in sites if k in tabs[O] and sel(r := tabs[O][k], mode) and r[7] in AA]
    obs_site = [persite[k] for k in sites if persite[k] in AA]
    for lab, obs in (('pooled', obs_pool), ('one-per-site', obs_site)):
        ll = {A: sum(lp(A, x) for x in obs) for A in list(AA) + ['bg']}
        rk = sorted(ll, key=ll.get, reverse=True)
        print(f'C. LL of {len(obs)} TAG-site control residues ({lab}); delta vs best: ' + '  '.join(f'{A}:{ll[A]-ll[rk[0]]:.1f}' for A in rk[:6])
              + f' | E-Q={ll["E"]-ll["Q"]:.1f} E-bg={ll["E"]-ll["bg"]:.1f} Q-bg={ll["Q"]-ll["bg"]:.1f}')
        if lab == 'one-per-site' and obs:
            bs = []
            for _ in range(2000):
                s = [random.choice(obs) for _ in obs]
                bs.append((sum(lp('E', x) - lp('Q', x) for x in s), sum(lp('E', x) - lp('bg', x) for x in s)))
            bs.sort(key=lambda t: t[0]); eq = (bs[50][0], bs[1949][0])
            bs.sort(key=lambda t: t[1]); eb = (bs[50][1], bs[1949][1])
            print(f'   bootstrap 95% CI (sites resampled): E-Q [{eq[0]:.1f}, {eq[1]:.1f}]  E-bg [{eb[0]:.1f}, {eb[1]:.1f}]')
    # D. mixture (EM) over amino-acid profiles + background, one residue per site
    comps = list(AA) + ['bg']; w = {A: 1 / len(comps) for A in comps}
    for _ in range(500):
        acc = collections.Counter()
        for x in obs_site:
            p = {A: w[A] * math.exp(lp(A, x)) for A in comps}; z = sum(p.values())
            for A in comps: acc[A] += p[A] / z
        w = {A: acc[A] / max(1, len(obs_site)) for A in comps}
    print('D. EM mixture weights (one residue per site): ' + ' '.join(f'{A}:{w[A]:.2f}' for A in sorted(w, key=w.get, reverse=True)[:6]))
