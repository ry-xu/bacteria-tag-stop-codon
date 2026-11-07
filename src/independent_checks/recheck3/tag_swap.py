"""Test-vs-test codon at aligned positions (from out/cols/<T>__<O>.tsv): at T's TAG columns, O's codon and residue;
calibration with T's GAG/GAA/CAG/CAA columns. Also lists every TAG site where O's codon differs from TAG.
Optional 4th arg: comma-separated controls -> extra set 'coding' = T columns that are also in the 'strict' set of a
T->control table with a sense (non-'*') control residue, i.e. the T frame is a protein-coding frame conserved at genus level
(A1 vs A2 are ~97% identical, so their alignments otherwise include non-coding / antisense frames).
Usage: python3 -I tag_swap.py <dir> T O [ctrl1,ctrl2,...]"""
import sys, collections
R, T, O = sys.argv[1:4]; ctrls = sys.argv[4].split(',') if len(sys.argv) > 4 else []
rows = [l.rstrip('\n').split('\t') for l in open(f'{R}/out/cols/{T}__{O}.tsv')][1:]
def sel(r, mode):
    if mode == 'census': return r[10] == '1'
    if mode == 'genome': return float(r[11]) >= 0.4 and int(r[12]) >= 10
    return float(r[11]) >= 0.5 and int(r[12]) >= 10 and int(r[13]) == 0
coding = set()
for c in ctrls:
    for l in list(open(f'{R}/out/cols/{T}__{c}.tsv'))[1:]:
        r = l.rstrip('\n').split('\t')
        if float(r[11]) >= 0.5 and int(r[12]) >= 10 and int(r[13]) == 0 and r[7] != '*': coding.add((r[0], r[1], r[2]))
_sel = sel
def sel(r, mode): return (r[0], r[1], r[2]) in coding and float(r[11]) >= 0.5 if mode == 'coding' else _sel(r, mode)
for mode in (['coding'] if ctrls else []) + ['census', 'strict', 'genome', 'all']:
    print(f'== {T} -> {O} [{mode}]')
    for cods in ('TAG', 'GAG', 'GAA', 'CAG', 'CAA'):
        cnt = collections.Counter(r[6] for r in rows if r[3] == cods and (mode == 'all' or sel(r, mode)))
        n = sum(cnt.values()); same = cnt[cods]
        oth = collections.Counter({k: v for k, v in cnt.items() if k != cods})
        print(f'   {cods}: n={n} identical={same} ({same/max(1,n):.0%}); changed -> ' + ' '.join(f'{k}:{v}' for k, v in oth.most_common(8)))
print('-- TAG columns where the other genome has a different codon (' + ('coding set' if ctrls else 'all') + '):')
for r in rows:
    if r[3] == 'TAG' and r[6] != 'TAG' and (not ctrls or sel(r, 'coding')): print('   ' + '\t'.join(r[i] for i in (0, 1, 2, 6, 7, 9, 10, 11, 12)))
