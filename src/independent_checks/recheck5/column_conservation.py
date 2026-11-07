"""Column conservation at T's TAG sites vs T's own codons, per control-consensus residue X in {W, Y, F}: fraction of outgroup
residues (loose alignment, J/S/E. coli/B. subtilis) equal to X, for T codon = TAG, = X's own codons, = TGG (W) and = anything
else. Sites: coding filter, >= 2 controls, unanimous (as tag_meaning5.py). Usage: python3 -I column_conservation.py <dir> T ctrls outgs"""
import sys, collections
sys.argv = sys.argv[:4] + [sys.argv[4], 'coding']
exec(open(sys.argv[0].replace('column_conservation.py', 'tag_meaning5.py')).read().split("print(f'##### {T}")[0])
keys = set().union(*[set(t) for t in tabs.values()]); cons = {}; tcod = {}
for k in keys:
    rs = [tabs[O][k][7] for O in ctrls if k in tabs[O]]
    if len(rs) >= 2 and len(set(rs)) == 1 and rs[0] in AA: cons[k] = rs[0]; tcod[k] = next(tabs[O][k][3] for O in ctrls if k in tabs[O])
ot = {O: load(O, strict=False) for O in outg}
own = {'W': {'TGG'}, 'Y': {'TAT', 'TAC'}, 'F': {'TTT', 'TTC'}}
for X in ('W', 'Y', 'F'):
    for lab, sel in (('TAG', lambda c: c == 'TAG'), (f'own({X})', lambda c: c in own[X]), ('TGG(W)', lambda c: c == 'TGG' and X != 'W'),
                     ('other', lambda c: c != 'TAG' and c not in own[X] and c != 'TGG')):
        ks = [k for k in cons if cons[k] == X and sel(tcod[k])]
        n = m = arom = 0; per = []
        for k in ks:
            rs = [ot[O][k][7] for O in outg if k in ot[O]]
            if not rs: continue
            n += len(rs); m += sum(r == X for r in rs); arom += sum(r in 'WYF' for r in rs); per.append(sum(r == X for r in rs) / len(rs))
        if ks: print(f'column {X}  T={lab:8s} sites={len(ks):5d} with_outgroup={len(per):5d}  outgroup residue == {X}: {m/max(1,n):.0%}  aromatic(WYF): {arom/max(1,n):.0%}  mean per-site {sum(per)/max(1,len(per)):.0%}')
# resampling test for W columns: mean per-site outgroup-W fraction of the TAG sites vs random draws of the same number of
# sites from T=TGG W-columns and from T=other W-columns
import random; random.seed(3)
def persite(sel):
    out = []
    for k in cons:
        if cons[k] != 'W' or not sel(tcod[k]): continue
        rs = [ot[O][k][7] for O in outg if k in ot[O]]
        if rs: out.append(sum(r == 'W' for r in rs) / len(rs))
    return out
tg, ow, oth = persite(lambda c: c == 'TAG'), persite(lambda c: c == 'TGG'), persite(lambda c: c not in ('TAG', 'TGG'))
m = sum(tg) / len(tg); N = 20000
p_oth = sum(sum(random.sample(oth, len(tg))) / len(tg) >= m for _ in range(N)) / N
p_tgg = sum(sum(random.sample(ow, len(tg))) / len(tg) <= m for _ in range(N)) / N
print(f'W columns: TAG per-site mean {m:.2f} (n={len(tg)}); P(random {len(tg)} T-substituted W-columns >= this) = {p_oth:.4f}; '
      f'P(random {len(tg)} TGG W-columns <= this) = {p_tgg:.3f}')
