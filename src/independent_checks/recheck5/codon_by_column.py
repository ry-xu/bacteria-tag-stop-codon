"""Codon-meaning contingency test. Same consensus sites as tag_meaning5.py (coding filter, >= 2 controls, unanimous residue).
For each control residue class X: counts of T codons TAG, TGG (W), CAR (Q), GAR (E), TAY (Y). If TAG encoded amino acid A,
TAG should be distributed over column classes like A's own codons (up to a constant). Prints the table, the per-class TAG
ratio to each candidate codon family, and multinomial log-likelihoods of the TAG column-class counts under 'TAG distributed
like TGG / CAR / GAR / TAY' (add-0.5), plus Fisher exact tests (2x2: TAG vs family at W columns vs Q columns).
Usage: python3 -I codon_by_column.py <dir> T ctrl1,ctrl2,..."""
import sys, collections, math, subprocess
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], '', 'coding']
src = open(sys.argv[0].replace('codon_by_column.py', 'tag_meaning5.py')).read()
src = src.split("print(f'##### {T}")[0]     # reuse loading/consensus code only
exec(src)
keys = set().union(*[set(t) for t in tabs.values()])
cons = {}; tcod = {}
for k in keys:
    rs = [tabs[O][k][7] for O in ctrls if k in tabs[O]]
    if len(rs) < 2 or len(set(rs)) != 1 or rs[0] not in AA: continue
    tcod[k] = next(tabs[O][k][3] for O in ctrls if k in tabs[O]); cons[k] = rs[0]
fam = {'TAG': {'TAG'}, 'TGG': {'TGG'}, 'CAR': {'CAG', 'CAA'}, 'GAR': {'GAG', 'GAA'}, 'TAY': {'TAT', 'TAC'}}
tab = {f: collections.Counter(cons[k] for k in cons if tcod[k] in cs) for f, cs in fam.items()}
print(f'{T}: unanimous coding sites {len(cons)}')
print('ctrl X  ' + ' '.join(f'{f:>6s}' for f in fam))
for X in sorted(AA, key=lambda X: -tab['TAG'][X] * 1e6 - tab['TGG'][X]):
    print(f'   {X}    ' + ' '.join(f'{tab[f][X]:6d}' for f in fam))
n = sum(tab['TAG'].values())
for f in ('TGG', 'CAR', 'GAR', 'TAY'):
    tot = sum(tab[f].values())
    ll = sum(tab['TAG'][X] * math.log((tab[f][X] + 0.5) / (tot + 10)) for X in AA)
    print(f'LL(TAG column classes | distributed like {f}) = {ll:.1f}')
def fisher(a, b, c, d):   # one-sided P(table at least as extreme in a's direction), hypergeometric
    n1, n2, k = a + b, c + d, a + c
    def p(x): return math.comb(n1, x) * math.comb(n2, k - x) / math.comb(n1 + n2, k)
    return sum(p(x) for x in range(a, min(n1, k) + 1))
for f in ('CAR', 'TGG'):
    a, b = tab['TAG']['W'], tab[f]['W']; c, d = tab['TAG']['Q'], tab[f]['Q']
    print(f'2x2 W-columns vs Q-columns: TAG {a} / {c}; {f} {b} / {d}; one-sided Fisher P(TAG enriched at W relative to {f}) = {fisher(a, b, c, d):.2g}')
