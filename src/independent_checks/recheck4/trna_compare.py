"""Compare CTA tRNAs of the test genomes with (a) each other, (b) the Christensenellales TAG-reader CTA tRNAs given in the brief,
(c) every Gln/Glu (and all other) tRNA in tests+controls: global identity (match/len of longer), plus mean identity per isotype
family; prints the structure lines (tRNAscan -f) for the CTA tRNAs and the Gln/Glu tRNAs of the same genome.
Usage: uv run --no-project --with biopython python -I trna_compare.py <dir> acc..."""
import sys, collections
from Bio import Align
R = sys.argv[1]; accs = sys.argv[2:]; D = f'{R}/out/trna'
CHR = {'Chris_reader_1': 'GGCGCGTTCGTCTAGTGGTCCAGGACAGCGGTCTCTAAAATCGCGAACATGGGTTCGAGCCCCATACGCGCTGCCA',
       'Chris_reader_2': 'CGCCCCATCCTCTAGCGGTCCAGGAGAGAGGTCTCTAAAACCTTGAACATCGGTTCGAATCCGATTGGGGTGACCA'}
CHR = {k: (v[:-3] if v.endswith('CCA') else v) for k, v in CHR.items()}
def rf(p):
    d = {}; k = None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): k = l[1:]; d[k] = ''
        else: d[k] += l.upper()
    return d
al = Align.PairwiseAligner(); al.mode = 'global'; al.match_score = 1; al.mismatch_score = -1; al.open_gap_score = -2; al.extend_gap_score = -1
def ident(x, y):
    a = al.align(x, y)[0]; m = sum(1 for (p, q) in zip(*a) if p == q and p != '-'); return m / max(len(x), len(y))
allt = {}
for a in accs:
    for k, v in rf(f'{D}/{a}.trna.fa').items():
        f = k.split(); allt[(a, f[0])] = (f[3] + f[4], v)
cta = {k: v for k, v in allt.items() if 'CTA' in v[0].upper()}
q = {f'{a}|{n}': v[1] for (a, n), v in cta.items()}; q.update(CHR)
names = list(q)
print('pairwise identity among CTA tRNAs (tests) and Christensenellales reader CTA tRNAs:')
for i in range(len(names)):
    for j in range(i + 1, len(names)): print(f'  {names[i]:45s} {names[j]:45s} {ident(q[names[i]], q[names[j]]):.3f}')
for n in names:
    fam = collections.defaultdict(list)
    for (a, nn), (iso, s) in allt.items():
        if f'{a}|{nn}' == n: continue
        fam[iso.split('(')[0]].append(ident(q[n], s))
    rk = sorted(fam, key=lambda k: -max(fam[k]))
    print(f'\n{n}: best identity per isotype family (max / mean over all tRNAs of that isotype in {len(accs)} genomes):')
    print('  ' + '  '.join(f'{k}:{max(fam[k]):.2f}/{sum(fam[k])/len(fam[k]):.2f}' for k in rk[:8]))
# structures
for a in accs:
    lines = open(f'{D}/{a}.struct.txt').read().split('\n\n')
    for blk in lines:
        if any(t in blk for t in ('Type: Sup', 'Type: Gln', 'Type: Glu')) and a in sys.argv[2:4]:
            print('\n' + a + '\n' + blk.strip())
