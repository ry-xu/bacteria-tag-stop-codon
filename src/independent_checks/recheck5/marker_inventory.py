"""Gene-family presence per genome from the census tblastn (best HSP e<1e-10, qcov>=0.5 per gene name), plus loci
(merged spans) per family. Prints families present in >=2 of the reference genomes but absent from T, and per-family
copy numbers that differ. Usage: python3 -I marker_inventory.py <dir> T ref1,ref2,..."""
import sys, collections
R, T, refs = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
def fam(a):
    L = collections.defaultdict(list)
    for l in open(f'{R}/out/{a}.tblastn.tsv'):
        f = l.rstrip('\n').split('\t'); q = f[0]
        if not (q.startswith('Ec|') or q.startswith('Bs|')): continue
        qs, qe, ss, se, ev, qlen = int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[8]), int(f[10])
        if ev >= 1e-10 or (qe - qs + 1) / qlen < 0.5: continue
        s = f[1].split('|')[1] if '|' in f[1] else f[1]
        L[q.split('|')[1]].append((s, min(ss, se), max(ss, se), float(f[2])))
    out = {}
    for g, hs in L.items():
        hs.sort(); m = []
        for s, a_, b, p in hs:
            if m and m[-1][0] == s and a_ <= m[-1][2]: m[-1][2] = max(m[-1][2], b)
            else: m.append([s, a_, b])
        out[g] = m
    return out
F = {a: fam(a) for a in [T] + refs}
allg = sorted(set().union(*[set(F[a]) for a in refs]))
miss = [g for g in allg if g not in F[T] and sum(g in F[a] for a in refs) >= 2]
print(f'{T}: families {len(F[T])}; refs: ' + ' '.join(f'{a}:{len(F[a])}' for a in refs))
print(f'absent from T but in >=2 refs ({len(miss)}): ' + ' '.join(miss))
extra = [g for g in F[T] if sum(g in F[a] for a in refs) == 0]
print(f'present in T, absent from all refs ({len(extra)}): ' + ' '.join(extra))
cn = [(g, len(F[T][g]), [len(F[a].get(g, [])) for a in refs]) for g in F[T] if len(F[T][g]) > max(len(F[a].get(g, [])) for a in refs)]
print('T copy number above every ref: ' + ' '.join(f'{g}:{n}{r}' for g, n, r in cn))
for g, n, r in cn:
    for s, a_, b in F[T][g]: print(f'   {g} {s}:{a_}-{b}')
