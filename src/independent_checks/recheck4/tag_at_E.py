"""List TAG sites of T whose controls unanimously show residue X (default E), with the sister's codon at that position,
the control codons and alignment quality. Usage: python3 -I tag_at_E.py <dir> T sister ctrl1,ctrl2,... [X]"""
import sys
R, T, S, ctrls = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(','); X = sys.argv[5] if len(sys.argv) > 5 else 'E'
C = f'{R}/out/cols'
def ok(f): return float(f[11]) >= 0.5 and int(f[12]) >= 10 and int(f[13]) == 0
def load(O):
    d = {}
    for l in list(open(f'{C}/{T}__{O}.tsv'))[1:]:
        f = l.rstrip('\n').split('\t')
        if ok(f): d[(f[0], f[1], f[2])] = f
    return d
tabs = {O: load(O) for O in ctrls}; sis = load(S)
for k in sorted(set().union(*tabs.values())):
    rs = [tabs[O][k] for O in ctrls if k in tabs[O]]
    if len(rs) < 2 or rs[0][3] != 'TAG' or len({r[7] for r in rs}) != 1 or rs[0][7] != X: continue
    s = sis.get(k)
    print(f'{k[0]}:{k[1]}:{k[2]}  ctrl codons ' + ','.join(r[6] for r in rs) + f'  win10_id ' + ','.join(r[11] for r in rs) +
          f'  sister: {s[6] + "(" + s[7] + ")" if s else "not aligned"}')
