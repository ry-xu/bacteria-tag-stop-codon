"""Collect tRNAscan-SE calls; print CTA-anticodon tRNAs, whether identical between genomes, Glu/Gln/Sup sets.
Writes out/trna/sup_cta.fa and out/trna/all_trna.fa. Usage: python3 -I trna_extract.py <dir> acc..."""
import sys
R = sys.argv[1] + '/out/trna'; accs = sys.argv[2:]
def rf(p):
    d = {}; k = None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): k = l[1:]; d[k] = ''
        else: d[k] += l
    return d
allt = {}
for a in accs:
    for k, v in rf(f'{R}/{a}.trna.fa').items(): allt[(a, k)] = v
sup = [(k, v) for k, v in allt.items() if '(CTA)' in k[1]]
for k, v in sup: print(k[0], k[1]); print(v)
print('distinct CTA tRNA sequences:', len({v for _, v in sup}))
open(f'{R}/sup_cta.fa', 'w').write(''.join(f'>{k[0]}|{k[1].split()[0]}\n{v}\n' for k, v in sup))
open(f'{R}/all_trna.fa', 'w').write(''.join(f'>{k[0]}|{k[1].split()[0]}|{k[1].split()[3]}{k[1].split()[4]}\n{v}\n' for k, v in allt.items()))
for a in accs:
    print(a, 'total', sum(1 for k in allt if k[0] == a), sorted(k[1].split()[3] + k[1].split()[4] for k in allt if k[0] == a and k[1].split()[3] in ('Glu', 'Gln', 'Sup')))
