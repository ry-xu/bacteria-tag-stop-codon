"""Summarise what TAG encodes, from out/cols/<T>__<O>.tsv tables (codon_columns.py).
Sets: 'census' = T codon inside a frame-matched conserved-protein HSP (the step-1 set);
      'genome' = all columns with win10_ident >= 0.4 and >= 10 columns from the HSP edge (genome-wide, more power).
A. TAG sites: residue in each same-genus control (counts), and per-site consensus over the controls.
B. Calibration: same statistic for test codons GAG/GAA (Glu), CAG/CAA (Gln) and other TAG neighbours.
C. Likelihood: profile P(control residue | test amino acid A) from all sense codons for A in the same set;
   log-likelihood of the TAG-site control residues under each A (pooled over controls; '*' columns excluded).
D. Test-vs-test codon swap: at TAG sites of one test genome, codon in the other test genome; calibration with GAG/GAA/CAG/CAA.
Usage: python3 -I tag_meaning.py <dir>"""
import sys, collections, math
R = sys.argv[1]; C = f'{R}/out/cols'
tests = ['GCA_024404275.1', 'GCA_045599485.1']; ctrls = ['GCA_934796175.1', 'GCA_004556155.1', 'GCA_963624325.1']
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
AA = 'ACDEFGHIKLMNPQRSTVWY'
def load(T, O):
    rows = [l.rstrip('\n').split('\t') for l in open(f'{C}/{T}__{O}.tsv')][1:]
    return {(r[0], r[1], r[2]): r for r in rows}
def sel(r, mode):
    return r[10] == '1' if mode == 'census' else (float(r[11]) >= 0.4 and int(r[12]) >= 10)
def fmt(cnt, k=8):
    n = sum(cnt.values()); return f'n={n} ' + ' '.join(f'{a}:{c}' for a, c in cnt.most_common(k))
for T in tests:
    tabs = {O: load(T, O) for O in ctrls}
    for mode in ('census', 'genome'):
        print(f'\n===== {T}  set={mode}')
        # A. TAG sites
        sites = set()
        for O in ctrls: sites |= {k for k, r in tabs[O].items() if r[3] == 'TAG' and sel(r, mode)}
        print(f'A. TAG sites with >=1 control aligned: {len(sites)}')
        for O in ctrls:
            cnt = collections.Counter(tabs[O][k][7] for k in sites if k in tabs[O] and sel(tabs[O][k], mode))
            print(f'   vs {O}: {fmt(cnt, 12)}')
        cons = collections.Counter()
        for k in sites:
            res = [tabs[O][k][7] for O in ctrls if k in tabs[O] and sel(tabs[O][k], mode)]
            top, n = collections.Counter(res).most_common(1)[0]
            cons[top if n > len(res) / 2 else 'mixed'] += 1
        print(f'   per-site majority over controls: {fmt(cons, 12)}')
        # B. calibration
        print('B. calibration (pooled over 3 controls; % = share of control residues):')
        for cods in (('GAG',), ('GAA',), ('CAG',), ('CAA',), ('GAT', 'GAC'), ('AAG', 'AAA'), ('TAT', 'TAC'), ('TTG',), ('TCG',), ('TGG',)):
            cnt = collections.Counter()
            for O in ctrls:
                for r in tabs[O].values():
                    if r[3] in cods and sel(r, mode): cnt[r[7]] += 1
            n = sum(cnt.values())
            print(f'   test {"/".join(cods)} ({CODE[cods[0]]}): n={n} ' + ' '.join(f'{a}:{cnt[a]/n:.1%}' for a, _ in cnt.most_common(6)))
        # TAG same format
        cnt = collections.Counter()
        for O in ctrls:
            for r in tabs[O].values():
                if r[3] == 'TAG' and sel(r, mode): cnt[r[7]] += 1
        n = sum(cnt.values())
        print(f'   test TAG          : n={n} ' + ' '.join(f'{a}:{cnt[a]/n:.1%}' for a, _ in cnt.most_common(8)))
        # C. likelihood
        prof = {A: collections.Counter() for A in AA}
        for O in ctrls:
            for r in tabs[O].values():
                a = CODE.get(r[3]);
                if a in prof and sel(r, mode) and r[7] in AA: prof[a][r[7]] += 1
        obs = collections.Counter()
        for O in ctrls:
            for r in tabs[O].values():
                if r[3] == 'TAG' and sel(r, mode) and r[7] in AA: obs[r[7]] += 1
        ll = {}
        for A in AA:
            n = sum(prof[A].values())
            ll[A] = sum(c * math.log((prof[A][x] + 0.5) / (n + 10)) for x, c in obs.items())
        rk = sorted(ll, key=ll.get, reverse=True)
        print(f'C. log-likelihood of {sum(obs.values())} TAG-site control residues under each amino-acid profile (best first, delta vs best):')
        print('   ' + '  '.join(f'{A}:{ll[A]-ll[rk[0]]:.1f}' for A in rk[:8]))
# D. test vs test
print('\n===== D. test-vs-test codon at aligned positions')
for T, O in ((tests[0], tests[1]), (tests[1], tests[0])):
    tab = load(T, O)
    for mode in ('census', 'genome'):
        out = []
        for cods in (('TAG',), ('GAG',), ('GAA',), ('CAG',), ('CAA',)):
            cnt = collections.Counter(r[6] for r in tab.values() if r[3] in cods and sel(r, mode))
            n = sum(cnt.values())
            same = cnt[cods[0]]
            other = collections.Counter({k: v for k, v in cnt.items() if k != cods[0]})
            out.append(f'   {T} {cods[0]} -> {O} codon [{mode}]: n={n}, identical={same}, changed: ' + ' '.join(f'{k}:{v}' for k, v in other.most_common(8)))
        print('\n'.join(out))
