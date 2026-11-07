"""recheck5 tRNA report: every tRNA-Trp(CCA), CTA-anticodon, Sup and Undet tRNA (tRNAscan-SE -B, --detail, -s) in the given
genomes, with sequence, score, isotype-model top scores; pairwise global identity (matches / longer length) among them;
anticodon-arm pairing of each sequence projected onto the T tRNA-Trp (5-bp anticodon stem, 7-nt loop, positions found by
aligning to T's tRNA-Trp, whose loop is known from tRNAscan); and a blastn search of T's tRNA-Trp against all genomes for any
further Trp-like copies tRNAscan did not report.
Usage: uv run --no-project --with biopython python -I trna_trp_cta.py <dir> T acc..."""
import sys, subprocess
from Bio import Align
R, T = sys.argv[1], sys.argv[2]; accs = sys.argv[3:]; D = f'{R}/out/trna'
def rf(p):
    d = {}; k = None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): k = l[1:]; d[k] = ''
        else: d[k] += l.upper()
    return d
def iso(a):
    out = {}; hdr = None
    for l in open(f'{D}/{a}.iso.txt'):
        f = [x.strip() for x in l.rstrip('\n').split('\t')]
        if f[0] == 'tRNAscanID': hdr = f; continue
        if hdr and len(f) == len(hdr): out[f[0]] = {h: float(v) for h, v in zip(hdr[2:], f[2:])}
    return out
al = Align.PairwiseAligner(); al.mode = 'global'; al.match_score = 1; al.mismatch_score = -1; al.open_gap_score = -2; al.extend_gap_score = -1
def ident(x, y):
    a = al.align(x, y)[0]; m = sum(1 for (p, q) in zip(*a) if p == q and p != '-'); return m / max(len(x), len(y))
sel = []
for a in [T] + accs:
    fa = rf(f'{D}/{a}.trna.fa'); I = iso(a)
    for k, v in fa.items():
        f = k.split(); name = f[0].split('::')[0]; typ, ac = f[3] if len(f) > 3 else '?', f[4] if len(f) > 4 else '?'
        if 'Trp' in typ or 'CTA' in ac.upper() or 'Sup' in typ or 'Undet' in typ:
            sc = I.get(name, {}); top = sorted(sc.items(), key=lambda x: -x[1])[:3]
            sel.append((a, name, typ + ac, k, v.replace('U', 'T'), top))
print('## Trp / CTA / Sup / Undet tRNAs')
for a, n, t, k, v, top in sel:
    print(f'{a} {n} {t} | {" ".join(k.split()[1:])}\n    seq {v}\n    top isotype models: ' + ' '.join(f'{x}={y:.1f}' for x, y in top))
# T's tRNA-Trp is the reference for numbering
tt = next(s for s in sel if s[0] == T and 'Trp' in s[2])[4]
# anticodon loop in T Trp: find CCA in the anticodon arm region (tRNAscan reported 33-35); loop 31-37, stem 26-30 / 38-42 (1-based)
ac = 33
print(f'\nreference T tRNA-Trp: anticodon arm {tt[25:30]}-[{tt[30:37]}]-{tt[37:42]} (stem 26-30, loop 31-37, anticodon {tt[32:35]})')
PAIR = {('G', 'C'), ('C', 'G'), ('A', 'T'), ('T', 'A'), ('G', 'T'), ('T', 'G')}
def project(seq):
    a = al.align(tt, seq)[0]; m = {}
    for (rs, re_), (qs, qe) in zip(*a.aligned):
        for k in range(re_ - rs): m[rs + k] = seq[qs + k]
    return m
print('\n## anticodon arm projected on T tRNA-Trp numbering (stem 26-30 : 38-42, loop 31-37); pairs listed outer->inner')
for a, n, t, k, v, top in sel:
    m = project(v)
    stem5 = ''.join(m.get(i, '-') for i in range(25, 30)); loop = ''.join(m.get(i, '-') for i in range(30, 37)); stem3 = ''.join(m.get(i, '-') for i in range(37, 42))
    pairs = [(stem5[i], stem3[4 - i]) for i in range(5)]
    print(f'{a} {n:28s} {t:10s} {stem5}-[{loop}]-{stem3}  anticodon {loop[2:5]}  pairs ' + ' '.join(f'{x}{y}{"" if (x, y) in PAIR else "!"}' for x, y in pairs))
print('\n## pairwise identity')
for i in range(len(sel)):
    for j in range(i + 1, len(sel)):
        print(f'{sel[i][0]}|{sel[i][1]}({sel[i][2]})\t{sel[j][0]}|{sel[j][1]}({sel[j][2]})\t{ident(sel[i][4], sel[j][4]):.3f}')
# blastn T tRNA-Trp vs all genomes
q = f'{D}/T_trp.fa'; open(q, 'w').write(f'>T_tRNA-Trp\n{tt}\n')
print('\n## blastn of T tRNA-Trp (word 7, e<1e-3) against each genome')
for a in [T] + accs:
    out = subprocess.run(['blastn', '-task', 'blastn-short', '-word_size', '7', '-query', q, '-db', f'{R}/data/{a}', '-evalue', '1e-3',
                          '-outfmt', '6 sseqid pident length qstart qend sstart send evalue'], capture_output=True, text=True, check=True).stdout
    for l in out.splitlines(): print(f'   {a} {l}')
