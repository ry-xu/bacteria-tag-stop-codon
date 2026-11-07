"""How divergent is T's RF1 compared with T's other conserved genes? For every census query (Ec|/Bs| conserved set plus the
RF queries Ec_prfA/Bs_prfA), take the best HSP (e<1e-10, qcov>=0.5) in T and in a reference genome, project both on query
coordinates and compute T-vs-reference identity over shared positions (>= 100). Reports the distribution (median, quartiles)
and the percentile rank of RF1 (prfA queries). Usage: python3 -I rf1_rate.py <dir> T ref..."""
import sys, statistics
R, T = sys.argv[1], sys.argv[2]; refs = sys.argv[3:]
def best(a):
    d = {}
    for l in open(f'{R}/out/{a}.tblastn.tsv'):
        f = l.rstrip('\n').split('\t'); q = f[0]
        if not (q.startswith('Ec|') or q.startswith('Bs|') or q.startswith('Ec_prfA') or q.startswith('Bs_prfA')): continue
        qs, qe, ev, bits, qlen = int(f[4]), int(f[5]), float(f[8]), float(f[9]), int(f[10])
        if ev >= 1e-10 or (qe - qs + 1) / qlen < 0.5: continue
        if q in d and d[q][0] >= bits: continue
        m = {}; p = qs
        for x, y in zip(f[13], f[14]):
            if x != '-':
                if y != '-': m[p] = y
                p += 1
        d[q] = (bits, m)
    return d
bt = best(T)
for r in refs:
    br = best(r); vals = {}
    for q in bt:
        if q not in br: continue
        sh = set(bt[q][1]) & set(br[q][1])
        if len(sh) >= 100: vals[q] = sum(bt[q][1][p] == br[q][1][p] for p in sh) / len(sh)
    v = sorted(vals.values()); rf = {q: x for q, x in vals.items() if 'prfA' in q}
    pct = {q: sum(1 for y in v if y < x) / len(v) for q, x in rf.items()}
    print(f'{T} vs {r}: genes={len(v)} identity median={statistics.median(v):.2f} IQR={v[len(v)//4]:.2f}-{v[3*len(v)//4]:.2f} | RF1: ' +
          ' '.join(f'{q.split("|")[0]}={x:.2f} (percentile {pct[q]:.0%})' for q, x in rf.items()))
