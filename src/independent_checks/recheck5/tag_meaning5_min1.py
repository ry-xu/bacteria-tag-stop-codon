"""Robustness variant of tag_meaning5.py: sites aligned (coding filter) in >= 1 control (residue = the controls' residue when
all aligned controls agree), profiles re-estimated on the same >= 1-control site set. LL per model, bootstrap CIs.
Usage: python3 -I tag_meaning5_min1.py <dir> T ctrl1,ctrl2,..."""
import sys, collections, math, random
sys.argv = sys.argv[:4] + ['', 'coding']
exec(open(sys.argv[0].replace('tag_meaning5_min1.py', 'tag_meaning5.py')).read().split("print(f'##### {T}")[0])
random.seed(9)
keys = set().union(*[set(t) for t in tabs.values()]); cons = {}; tcod = {}
for k in keys:
    rs = [tabs[O][k][7] for O in ctrls if k in tabs[O]]
    if len(set(rs)) == 1 and rs[0] in AA: cons[k] = rs[0]; tcod[k] = next(tabs[O][k][3] for O in ctrls if k in tabs[O])
prof = {A: collections.Counter() for A in AA}; bg = collections.Counter()
for k in cons:
    a = CODE.get(tcod[k])
    if a in prof: prof[a][cons[k]] += 1; bg[cons[k]] += 1
prof['bg'] = bg
def lp(A, x): n = sum(prof[A].values()); return math.log((prof[A][x] + 0.5) / (n + 10))
obs = [cons[k] for k in cons if tcod[k] == 'TAG']
print(f'{T}: >=1-control agreeing sites {len(cons)}; TAG sites {len(obs)}: ' + ' '.join(f'{a}:{n}' for a, n in collections.Counter(obs).most_common()))
M = list(AA) + ['bg']; ll = {A: sum(lp(A, x) for x in obs) for A in M}; rk = sorted(ll, key=ll.get, reverse=True)
print('LL delta vs best: ' + ' '.join(f'{A}:{ll[A]-ll[rk[0]]:.1f}' for A in rk[:5]) + f' | W-Q={ll["W"]-ll["Q"]:.1f} W-E={ll["W"]-ll["E"]:.1f} W-Y={ll["W"]-ll["Y"]:.1f} W-bg={ll["W"]-ll["bg"]:.1f}')
B = []
for _ in range(2000):
    s = [random.choice(obs) for _ in obs]; l = {A: sum(lp(A, x) for x in s) for A in M}
    B.append((l['W'] - l['Q'], l['W'] - l['Y'], l['W'] - l['bg'], max(l, key=l.get)))
for j, lab in enumerate(('W-Q', 'W-Y', 'W-bg')):
    v = sorted(x[j] for x in B); print(f'  bootstrap 95% CI {lab}: [{v[50]:.1f}, {v[1949]:.1f}]')
print('  bootstrap best: ' + ' '.join(f'{a}:{n/20:.1f}%' for a, n in collections.Counter(x[3] for x in B).most_common(3)))
