"""tRNA report (tRNAscan-SE -B --detail, -s isotype scores). For each genome: CTA-anticodon tRNAs (coords, score, note,
isotype-model scores and top isotype), sequence; Glu/Gln/Sup tRNA inventory. Then: identity of each CTA tRNA to every
other tRNA in all genomes (top hits -> likely source isotype), pairwise identity between CTA tRNAs of different genomes.
Usage: uv run --no-project --with biopython python -I trna_report.py <dir> acc..."""
import sys, re, collections
from Bio import Align
R = sys.argv[1]; accs = sys.argv[2:]; D = f'{R}/out/trna'
def rf(p):
    d = {}; k = None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): k = l[1:]; d[k] = ''
        else: d[k] += l.upper()
    return d
allt = {}; cta = []
for a in accs:
    fa = rf(f'{D}/{a}.trna.fa')
    iso = {}
    lines = [l.rstrip('\n') for l in open(f'{D}/{a}.iso.txt')]
    hdr = None
    for l in lines:
        f = l.split('\t')
        if hdr is None and len(f) > 3 and f[0].strip().lower().startswith('tRNAscanID'.lower()): hdr = [x.strip() for x in f]; continue
        if hdr and len(f) == len(hdr): iso[f[0].strip()] = dict(zip(hdr, [x.strip() for x in f]))
    for k, v in fa.items():
        name = k.split()[0]; m = re.search(r'(\w+)-(\w+)\s+\(([^)]+)\)\s+(\d+)\s+bp\s+Sc:\s+([\d.]+)', k)
        aa, ac = m.group(2) if False else k.split()[3], k.split()[4]
        allt[(a, name)] = (v, k)
        if '(cta)' in k.lower():
            cta.append((a, name, k, v))
    inv = collections.Counter(f'{k[1].split()[3]}{k[1].split()[4]}' for (aa_, n), k in allt.items() if aa_ == a)
    print(f'## {a}: {sum(inv.values())} tRNAs; Glu/Gln/Sup/Undet: ' + ' '.join(f'{x}:{c}' for x, c in sorted(inv.items()) if re.match(r'(Glu|Gln|Sup|Undet|SeC|Pyl)', x)))
    for (aa_, n), (v, k) in allt.items():
        if aa_ == a and '(cta)' in k.lower():
            print('   CTA tRNA:', k); print('   seq:', v)
            r = iso.get(n.split('::')[0]) or iso.get(n)
            if r:
                sc = {kk: vv for kk, vv in r.items() if kk not in ('tRNAscanID', 'anticodon', 'Anticodon')}
                print('   isotype row:', '  '.join(f'{kk}={vv}' for kk, vv in sc.items()))
print(f'\nCTA tRNAs found: {len(cta)}; distinct sequences: {len({c[3] for c in cta})}')
al = Align.PairwiseAligner(); al.mode = 'global'; al.match_score = 1; al.mismatch_score = -1; al.open_gap_score = -2; al.extend_gap_score = -1
def ident(x, y):
    a = al.align(x, y)[0]; m = sum(1 for (p, q) in zip(*a) if p == q and p != '-'); return m / max(len(x), len(y))
for i in range(len(cta)):
    for j in range(i + 1, len(cta)):
        print(f'CTA vs CTA: {cta[i][0]}|{cta[i][1]} vs {cta[j][0]}|{cta[j][1]} identity {ident(cta[i][3], cta[j][3]):.3f}')
for a, n, k, v in cta:
    res = sorted(((ident(v, s), aa_, kk.split()[3] + kk.split()[4], nn) for (aa_, nn), (s, kk) in allt.items() if (aa_, nn) != (a, n)), reverse=True)
    print(f'\nTop matches of {a}|{n} among all tRNAs:')
    for r in res[:10]: print(f'   {r[0]:.3f} {r[1]} {r[3]} {r[2]}')
