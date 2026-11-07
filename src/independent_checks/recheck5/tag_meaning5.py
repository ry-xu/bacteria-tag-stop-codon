"""recheck5: what amino acid does TAG encode in test genome T? (adapted from recheck4 tag_meaning4.py; same column tables)
Input: codon-column tables out/cols/<T>__<O>.tsv made by codon_columns.py (T ORFs with only TAA/TGA as stops, TAG -> X so the
aligned residue is not pulled toward any amino acid; tblastn vs O's raw DNA; frame competition; ungapped +-3).
Column filter ('strict'): win10_ident >= 0.5, >= 10 columns from the HSP edge, no stop in O's translated HSP within +-20 cols.

 1. Consensus: sites of T aligned (strict) in >= 2 controls; 'unanimous' if all controls show the same residue X.
    Table per T codon class -> distribution of X (calibration rows for sense codons, then TAG).
 2. Likelihood: P(X | T amino acid A) estimated from all T sense codons of A at unanimous sites (add-0.5 smoothing);
    'bg' = all sense sites pooled; 'mutnull' = TAG arose by a single-nucleotide change from a random sense codon at a site
    that is not under selection for the TAG meaning: mixture over T sense codons c one step from TAG, weighted by
    count(c) x relative rate (transitions kappa x transversions), of P(X | aa(c)). LL of the TAG-site residues under every A,
    bg and mutnull; bootstrap (2000x, resampling TAG sites) 95% CIs for LL(W)-LL(Q), LL(W)-LL(E), LL(W)-LL(bg),
    LL(W)-LL(mutnull), LL(W)-LL(best other).
 3. Mixture of W, Q, E, bg profiles: ML weights and 95% profile-likelihood interval for the W weight.
 4. Route check: TAG is one transition from TGG (W) and CAG (Q). Under 'TAG = Gln' TAGs should accumulate at Gln columns via
    CAG->TAG (common codon); under 'TAG = Trp' only at Trp columns. Reports the observed TAG count at unanimous-W vs unanimous-Q
    columns, and the expectation under each hypothesis from the numbers of TGG-at-W and CAG/CAA-at-Q sites.
 5. Conservation of the Trp columns (only if outgroup tables given): for T sites at control-unanimous W columns, the residue
    in each outgroup genome (more distant) at the same T codon; fraction of outgroups that also show W. Compared for
    T codon = TAG, T codon = TGG, and T codon != TGG/TAG (W columns where T itself tolerates a substitution).
Filter modes (env-free, 5th argument): 'strict' (above) or 'coding' = strict AND the control codon lies in-frame inside a
pyrodigal gene of that control (out/genes/<O>.genes.tsv from ctrl_genes.py; removes antisense/shadow-frame and non-coding
alignments, and every column where the control itself has a stop codon).
Usage: python3 -I tag_meaning5.py <dir> T ctrl1,ctrl2,... [outgroup1,outgroup2,... or ''] [strict|coding]"""
import sys, collections, math, random, bisect
R, T, ctrls = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
outg = sys.argv[4].split(',') if len(sys.argv) > 4 and sys.argv[4] else []
MODE = sys.argv[5] if len(sys.argv) > 5 else 'strict'
COMPL = str.maketrans('ACGT', 'TGCA')
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
def coding_checker(O):
    genes = {}
    for l in open(f'{R}/out/genes/{O}.genes.tsv'):
        c, b_, e_, st = l.split('\t'); genes.setdefault(c, []).append((int(b_), int(e_), int(st)))
    for c in genes: genes[c].sort()
    starts = {c: [g[0] for g in v] for c, v in genes.items()}
    gO = read_contigs(O)
    def in_gene(c, a, z, strand):
        v = genes.get(c, ()); i = bisect.bisect_right(starts.get(c, []), a)
        for b_, e_, st in v[max(0, i - 40):i]:
            if st != strand or a < b_ or z > e_: continue
            if (strand == 1 and (a - b_) % 3 == 0) or (strand == -1 and (e_ - z) % 3 == 0): return True
        return False
    def chk(f):
        oc, op, ocod = f[4], int(f[5]), f[6]; s_ = gO[oc]
        if s_[op-1:op+2] == ocod and in_gene(oc, op, op + 2, 1): return True
        if op >= 3 and s_[op-3:op].translate(COMPL)[::-1] == ocod and in_gene(oc, op - 2, op, -1): return True
        return False
    return chk
C = f'{R}/out/cols'; random.seed(5)
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
AA = 'ACDEFGHIKLMNPQRSTVWY'
def ok(f): return float(f[11]) >= 0.5 and int(f[12]) >= 10 and int(f[13]) == 0
def load(O, strict=True):
    d = {}; chk = coding_checker(O) if (strict and MODE == 'coding') else None
    for l in list(open(f'{C}/{T}__{O}.tsv'))[1:]:
        f = l.rstrip('\n').split('\t')
        if not strict or (ok(f) and (chk is None or chk(f))): d[(f[0], f[1], f[2])] = f
    return d
tabs = {O: load(O) for O in ctrls}
keys = set().union(*[set(t) for t in tabs.values()])
cons = {}; tcod = {}
for k in keys:
    rs = [tabs[O][k][7] for O in ctrls if k in tabs[O]]
    if len(rs) < 2: continue
    tcod[k] = next(tabs[O][k][3] for O in ctrls if k in tabs[O])
    cons[k] = rs[0] if len(set(rs)) == 1 else '~'
print(f'##### {T}  controls={",".join(ctrls)}  outgroups={",".join(outg) or "-"}  filter={MODE}')
ntag_any = sum(1 for k in keys if next(tabs[O][k][3] for O in ctrls if k in tabs[O]) == 'TAG')
print(f'1. strict-aligned T sites (any control): {len(keys)}, TAG among them: {ntag_any}; aligned in >=2 controls: {len(cons)} '
      f'(TAG {sum(1 for k in cons if tcod[k] == "TAG")}); unanimous: {sum(1 for v in cons.values() if v != "~")} '
      f'(TAG {sum(1 for k in cons if tcod[k] == "TAG" and cons[k] != "~")})')
def row(lab, sel):
    ks = [k for k in cons if sel(tcod[k])]; n = len(ks); u = [cons[k] for k in ks if cons[k] != '~']; cu = collections.Counter(u)
    print(f'   T {lab:8s} sites={n:6d} unanimous={len(u):6d} ({len(u)/max(1,n):.0%}) | ' + ' '.join(f'{a}:{c/max(1,len(u)):.1%}' for a, c in cu.most_common(6)))
for lab, cods in (('TGG(W)', {'TGG'}), ('CAG', {'CAG'}), ('CAA', {'CAA'}), ('GAG', {'GAG'}), ('GAA', {'GAA'}), ('TAY(Y)', {'TAT', 'TAC'}),
                  ('TTR(L)', {'TTA', 'TTG'}), ('AAR(K)', {'AAG', 'AAA'}), ('TCG(S)', {'TCG'}), ('TGY(C)', {'TGT', 'TGC'}), ('TAG', {'TAG'})):
    row(lab, lambda c, cods=cods: c in cods)
uni = [k for k in cons if cons[k] in AA]
prof = {A: collections.Counter() for A in AA}; bg = collections.Counter(); ccount = collections.Counter()
for k in uni:
    a = CODE.get(tcod[k])
    if a in prof: prof[a][cons[k]] += 1; bg[cons[k]] += 1; ccount[tcod[k]] += 1
prof['bg'] = bg
def lp(A, x): n = sum(prof[A].values()); return math.log((prof[A][x] + 0.5) / (n + 10))
# mutational null: one-step neighbours of TAG among sense codons
def ts(a, b): return {a, b} in ({'A', 'G'}, {'C', 'T'})
KAPPA = 3.0
nb = {}
for c, n in ccount.items():
    d = [(p, c[p], 'TAG'[p]) for p in range(3) if c[p] != 'TAG'[p]]
    if len(d) == 1: nb[c] = n * (KAPPA if ts(d[0][1], d[0][2]) else 1.0)
Z = sum(nb.values())
mut = {x: sum(w / Z * math.exp(lp(CODE[c], x)) for c, w in nb.items()) for x in AA}
print(f'   mutational-null neighbour weights (kappa={KAPPA}): ' + ' '.join(f'{c}({CODE[c]}):{w/Z:.2f}' for c, w in sorted(nb.items(), key=lambda t: -t[1])))
obs = [cons[k] for k in uni if tcod[k] == 'TAG']
tc = collections.Counter(obs)
print(f'   TAG unanimous residues (n={len(obs)}): ' + ' '.join(f'{a}:{n}' for a, n in tc.most_common()))
if not obs: sys.exit()
def LL(A, xs): return sum(math.log(mut[x]) if A == 'mutnull' else lp(A, x) for x in xs)
models = list(AA) + ['bg', 'mutnull']
ll = {A: LL(A, obs) for A in models}
rk = sorted(ll, key=ll.get, reverse=True)
print('2. LL of TAG-site residues, delta vs best: ' + '  '.join(f'{A}:{ll[A]-ll[rk[0]]:.1f}' for A in rk[:6]) +
      f' | W-Q={ll["W"]-ll["Q"]:.1f} W-E={ll["W"]-ll["E"]:.1f} W-bg={ll["W"]-ll["bg"]:.1f} W-mutnull={ll["W"]-ll["mutnull"]:.1f}')
B = []
for _ in range(2000):
    s = [random.choice(obs) for _ in obs]; l = {A: LL(A, s) for A in models}
    B.append((l['W'] - l['Q'], l['W'] - l['E'], l['W'] - l['bg'], l['W'] - l['mutnull'], l['W'] - max(v for A, v in l.items() if A != 'W'),
              max(l, key=l.get)))
for j, lab in enumerate(('W-Q', 'W-E', 'W-bg', 'W-mutnull', 'W-best_other')):
    v = sorted(x[j] for x in B); print(f'   bootstrap 95% CI {lab}: [{v[50]:.1f}, {v[1949]:.1f}]')
print('   bootstrap best model: ' + ' '.join(f'{a}:{n/20:.1f}%' for a, n in collections.Counter(x[5] for x in B).most_common(4)))
# leave-one-out influence: does a single site decide it?
loo = sorted((LL('W', obs[:i] + obs[i+1:]) - max(LL(A, obs[:i] + obs[i+1:]) for A in models if A != 'W'), obs[i]) for i in range(len(obs)))
print(f'   leave-one-out W-best_other range: [{loo[0][0]:.1f} (dropping {loo[0][1]}), {loo[-1][0]:.1f} (dropping {loo[-1][1]})]')
# 3. mixture W / Q / E / bg
P = {A: {x: math.exp(lp(A, x)) for x in AA} for A in ('W', 'Q', 'E', 'bg')}
grid = [i / 50 for i in range(51)]; best = None; byw = {}
for w in grid:
    for q in grid:
        for e in grid:
            if w + q + e > 1 + 1e-9: continue
            g_ = 1 - w - q - e
            L = sum(math.log(w * P['W'][x] + q * P['Q'][x] + e * P['E'][x] + g_ * P['bg'][x] + 1e-300) for x in obs)
            if w not in byw or L > byw[w][0]: byw[w] = (L, q, e)
            if best is None or L > best[0]: best = (L, w, q, e)
ci = [w for w in grid if byw[w][0] >= best[0] - 1.92]
print(f'3. mixture ML weights: W={best[1]:.2f} Q={best[2]:.2f} E={best[3]:.2f} bg={1-best[1]-best[2]-best[3]:.2f}; 95% profile interval for W weight [{min(ci):.2f}, {max(ci):.2f}]')
# 4. route check
nW_tgg = sum(1 for k in uni if cons[k] == 'W' and tcod[k] == 'TGG'); nQ_cag = sum(1 for k in uni if cons[k] == 'Q' and tcod[k] == 'CAG')
nQ_caa = sum(1 for k in uni if cons[k] == 'Q' and tcod[k] == 'CAA')
print(f'4. route check: unanimous-W columns with T=TGG: {nW_tgg}, with T=TAG: {tc["W"]}; unanimous-Q columns with T=CAG: {nQ_cag}, CAA: {nQ_caa}, T=TAG: {tc["Q"]}')
# audit fix: conditional binomial instead of a Poisson with a plug-in expectation. Given the TAG sites at W or Q columns,
# under 'TAG is used like CAG at Q columns as like TGG at W columns' the share at Q columns would be nQ_cag/(nQ_cag+nW_tgg).
nWQ = tc['W'] + tc['Q']; pq = nQ_cag / (nQ_cag + nW_tgg)
pbin = sum(math.comb(nWQ, i) * pq ** i * (1 - pq) ** (nWQ - i) for i in range(tc['Q'] + 1))
print(f'   TAG at W+Q columns = {nWQ}; expected share at Q columns if TAG were used like CAG (Gln) as it is like TGG (Trp): {pq:.2f}; '
      f'observed {tc["Q"]}/{nWQ}; conditional binomial P(<= obs) = {pbin:.2g} (assumes equal relative use, see text)')
# 5. conservation of W columns in outgroups
if outg:
    ot = {O: load(O, strict=False) for O in outg}
    print(f'5. W-column conservation in outgroups {",".join(outg)} (residue at the same T codon; loose: any aligned column)')
    for lab, sel in (('TAG', lambda c: c == 'TAG'), ('TGG', lambda c: c == 'TGG'), ('other(T substituted W)', lambda c: c not in ('TAG', 'TGG'))):
        ks = [k for k in uni if cons[k] == 'W' and sel(tcod[k])]
        nal = 0; nW = 0; allW = 0; anyal = 0; res = collections.Counter()
        for k in ks:
            rs = [ot[O][k][7] for O in outg if k in ot[O]]
            if not rs: continue
            anyal += 1; nal += len(rs); nW += sum(r_ == 'W' for r_ in rs); allW += all(r_ == 'W' for r_ in rs); res.update(rs)
        print(f'   T {lab:24s} W-columns={len(ks):5d} with outgroup data={anyal:5d}  outgroup residues W={nW/max(1,nal):.1%}  '
              f'all-outgroups-W={allW/max(1,anyal):.1%} | ' + ' '.join(f'{a}:{n}' for a, n in res.most_common(5)))
    # per-site list for TAG at W
    print('   TAG sites at control-unanimous W: site, control codons, outgroup residues')
    for k in sorted(k for k in uni if cons[k] == 'W' and tcod[k] == 'TAG'):
        print(f'     {k[0]}:{k[1]}:{k[2]}  ctrl ' + ','.join(tabs[O][k][6] for O in ctrls if k in tabs[O]) + '  outg ' +
              ' '.join(f'{O[4:13]}={ot[O][k][7] if k in ot[O] else "."}' for O in outg))
print('   all TAG sites with unanimous residue: ' + ' '.join(f'{k[0][-6:]}:{k[1]}={cons[k]}' for k in sorted(k for k in uni if tcod[k] == 'TAG')))
