"""recheck4: what does TAG encode? Independent analyses on the codon-column tables made by codon_columns.py
(test ORFs read with only TAA/TGA as stops and TAG -> X (neutral), tblastn vs other genome's raw DNA; one column per test codon;
frame-competition and ungapped +-3 filters). Column filter used throughout ('strict'): win10 identity >= 0.5, >= 10 columns from
the HSP edge, no stop in the other genome's translated HSP within +-20 columns.
 1. Consensus columns: test codon sites aligned in >= 2 controls that all show the SAME residue X. Table: test codon class ->
    distribution of X. Calibration rows for sense codons; TAG row.
 2. Profile matching (reverse view): for every unanimous residue X, enrichment of each test codon class at X-columns relative
    to all consensus columns; Pearson r between the TAG enrichment vector (log2, 20 aa) and each amino acid's (codon-family)
    vector -> which amino acid's codons TAG behaves like.
 3. Likelihood: P(consensus residue | test amino acid A) from all sense codons of A (add-0.5), background = all sense
    columns; LL of TAG-site consensus residues under Q, E, bg, every A; bootstrap over sites.
 4. Sister comparison (no controls): sister genome's codon at each TAG site of T (T__S table), vs the same for CAG/CAA/GAG/GAA;
    shared TAG sites (TAG in both at the orthologous position) vs the convergence expectation.
Usage: python3 -I tag_meaning4.py <dir> T sister ctrl1,ctrl2,..."""
import sys, collections, math, random
R, T, S, ctrls = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(','); C = f'{R}/out/cols'
random.seed(11)
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
AA = 'ACDEFGHIKLMNPQRSTVWY'
def ok(f): return float(f[11]) >= 0.5 and int(f[12]) >= 10 and int(f[13]) == 0
def load(T, O):
    d = {}
    for l in list(open(f'{C}/{T}__{O}.tsv'))[1:]:
        f = l.rstrip('\n').split('\t')
        if ok(f): d[(f[0], f[1], f[2])] = f
    return d
tabs = {O: load(T, O) for O in ctrls}
keys = set().union(*[set(t) for t in tabs.values()])
cons = {}; tcod = {}
for k in keys:
    rs = [tabs[O][k][7] for O in ctrls if k in tabs[O]]
    if len(rs) < 2: continue
    tcod[k] = next(tabs[O][k][3] for O in ctrls if k in tabs[O])
    cons[k] = rs[0] if len(set(rs)) == 1 else '~'     # '~' = controls disagree
cls = lambda c: 'TAG' if c == 'TAG' else CODE.get(c, '?')
print(f'##### {T}  controls={",".join(ctrls)}  sister={S}')
print(f'1. sites aligned in >=2 controls: {len(cons)}; unanimous: {sum(1 for v in cons.values() if v != "~")}')
def row(lab, sel):
    ks = [k for k in cons if sel(tcod[k])]; n = len(ks); u = [cons[k] for k in ks if cons[k] != '~']; cu = collections.Counter(u)
    print(f'   test {lab:8s} sites={n:6d} unanimous={len(u):6d} ({len(u)/max(1,n):.0%}) | among unanimous: ' + ' '.join(f'{a}:{c/max(1,len(u)):.1%}' for a, c in cu.most_common(6)))
    return cu
for lab, cods in (('CAG', {'CAG'}), ('CAA', {'CAA'}), ('GAG', {'GAG'}), ('GAA', {'GAA'}), ('AAR(K)', {'AAG', 'AAA'}), ('GAY(D)', {'GAT', 'GAC'}),
                  ('CGN+AGR', {'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'}), ('TGG(W)', {'TGG'}), ('TAY(Y)', {'TAT', 'TAC'}), ('TTR(L)', {'TTA', 'TTG'}), ('TAG', {'TAG'})):
    cu = row(lab, lambda c, cods=cods: c in cods)
tag_cu = collections.Counter(cons[k] for k in cons if tcod[k] == 'TAG' and cons[k] != '~')
# 2. profile matching
uni = [k for k in cons if cons[k] != '~' and cons[k] in AA]
nX = collections.Counter(cons[k] for k in uni); N = len(uni)
def enr(sel):
    ks = [k for k in uni if sel(tcod[k])]; n = len(ks); c = collections.Counter(cons[k] for k in ks)
    return {X: math.log2(((c[X] + 0.5) / (n + 10)) / (nX[X] / N)) for X in AA}, n
vt, nt = enr(lambda c: c == 'TAG')
def pear(a, b):
    xs = [a[X] for X in AA]; ys = [b[X] for X in AA]; mx = sum(xs) / 20; my = sum(ys) / 20
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
rs = {}
for A in AA:
    vA, nA = enr(lambda c, A=A: CODE.get(c) == A); rs[A] = (pear(vt, vA), nA)
print(f'2. TAG enrichment (log2 obs/exp) by unanimous control residue (n_TAG={nt}): ' + ' '.join(f'{X}:{vt[X]:+.1f}' for X in sorted(AA, key=lambda X: -vt[X])[:6]))
print('   Pearson r of TAG enrichment vector vs each amino acid codon family: ' + ' '.join(f'{A}:{rs[A][0]:.2f}' for A in sorted(AA, key=lambda A: -rs[A][0])[:6])
      + f' ... E:{rs["E"][0]:.2f}  (Q-family n={rs["Q"][1]}, E-family n={rs["E"][1]})')
# 3. likelihood on unanimous residues
prof = {A: collections.Counter() for A in AA}; bg = collections.Counter()
for k in uni:
    a = CODE.get(tcod[k])
    if a in prof: prof[a][cons[k]] += 1; bg[cons[k]] += 1
prof['bg'] = bg
def lp(A, x): n = sum(prof[A].values()); return math.log((prof[A][x] + 0.5) / (n + 10))
obs = [cons[k] for k in uni if tcod[k] == 'TAG']
ll = {A: sum(lp(A, x) for x in obs) for A in list(AA) + ['bg']}
rk = sorted(ll, key=ll.get, reverse=True)
print(f'3. LL of {len(obs)} unanimous TAG-site residues; delta vs best: ' + '  '.join(f'{A}:{ll[A]-ll[rk[0]]:.1f}' for A in rk[:5])
      + f' | LL(Q)-LL(E)={ll["Q"]-ll["E"]:.1f}  LL(Q)-LL(bg)={ll["Q"]-ll["bg"]:.1f}  LL(E)-LL(bg)={ll["E"]-ll["bg"]:.1f}')
if obs:
    bs = sorted(sum(lp('Q', x) - lp('E', x) for x in [random.choice(obs) for _ in obs]) for _ in range(2000))
    print(f'   bootstrap 95% CI LL(Q)-LL(E): [{bs[50]:.1f}, {bs[1949]:.1f}]')
    # expected unanimous-Q share if TAG were Q / E / bg
    q = tag_cu['Q'] / max(1, sum(tag_cu.values()))
    print(f'   unanimous-Q share at TAG sites {q:.1%} vs expected if Gln {prof["Q"]["Q"]/sum(prof["Q"].values()):.1%}, if Glu {prof["E"]["Q"]/sum(prof["E"].values()):.1%}, '
          f'background {bg["Q"]/sum(bg.values()):.1%}; unanimous-E share {tag_cu["E"]/max(1,sum(tag_cu.values())):.1%} vs if Gln {prof["Q"]["E"]/sum(prof["Q"].values()):.1%}, if Glu {prof["E"]["E"]/sum(prof["E"].values()):.1%}')
# 4. sister codons
st = load(T, S)
def sis(sel):
    c = collections.Counter(f[6] for k, f in st.items() if sel(f[3])); n = sum(c.values()); return c, n
print(f'4. sister {S}: codon at the orthologous position (strict columns)')
for lab, cods in (('TAG', {'TAG'}), ('CAG', {'CAG'}), ('CAA', {'CAA'}), ('GAG', {'GAG'}), ('GAA', {'GAA'}), ('all', None)):
    c, n = sis(lambda x, cods=cods: cods is None or x in cods)
    grp = lambda xs: sum(c[x] for x in xs) / max(1, n)
    print(f'   T {lab:4s} n={n:6d}: sister TAG {grp(["TAG"]):.1%}  CAR(Q) {grp(["CAG","CAA"]):.1%}  GAR(E) {grp(["GAG","GAA"]):.1%}  AAR(K) {grp(["AAG","AAA"]):.1%}  '
          f'other {1-grp(["TAG","CAG","CAA","GAG","GAA","AAG","AAA"]):.1%} | top: ' + ' '.join(f'{x}:{m}' for x, m in c.most_common(5)))
# 5. three-component mixture (Gln profile, Glu profile, background) on unanimous TAG-site residues: ML weights and
#    profile-likelihood 95% interval for the Glu weight.
if obs:
    pq = {x: math.exp(lp('Q', x)) for x in AA}; pe = {x: math.exp(lp('E', x)) for x in AA}; pb = {x: math.exp(lp('bg', x)) for x in AA}
    best = None; grid = [i / 100 for i in range(101)]
    tab = {}
    for we in grid:
        bl = None
        for wb in grid:
            if we + wb > 1: break
            L = sum(math.log((1 - we - wb) * pq[x] + we * pe[x] + wb * pb[x] + 1e-300) for x in obs)
            if bl is None or L > bl[0]: bl = (L, wb)
        tab[we] = bl
        if best is None or bl[0] > best[0]: best = (bl[0], we, bl[1])
    ci = [we for we in grid if tab[we][0] >= best[0] - 1.92]
    print(f'5. mixture on unanimous TAG residues: ML weights Gln={1-best[1]-best[2]:.2f} Glu={best[1]:.2f} background={best[2]:.2f}; '
          f'95% profile-likelihood interval for the Glu weight: [{min(ci):.2f}, {max(ci):.2f}]; LL(pure Glu+bg best) - LL(ML) = {max(tab[we][0] for we in grid if we >= 0.5) - best[0]:.1f} (Glu weight >= 0.5)')
