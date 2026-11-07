"""Release-factor comparison. Proteins: every std-mode RF1/RF2 HMM call >= 300 aa from out/rf_hmm.tsv (rf_candidates.faa) +
E. coli / B. subtilis prfA/prfB (data/rf_queries.faa). For each: DIAMOND blastp top-5 Swiss-Prot hits (reciprocal class check:
are the best hits RF1 (prfA) or RF2 (prfB)?), identity to the best Swiss-Prot RF1 and best RF2; pairwise local identities
(BLOSUM62, gap -11/-1; identity over aligned columns and coverage of the shorter protein); motif windows mapped by alignment
to E. coli RF1/RF2 (PxT/SPF loop +-6, GGQ +-4) printed side by side.
Usage: uv run --no-project --with biopython python -I rf_compare.py <dir>"""
import sys, subprocess, os
from Bio import Align
from Bio.Align import substitution_matrices
R = sys.argv[1]
def read_fasta(p):
    d = {}; k = None; buf = []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if k: d[k] = ''.join(buf)
            k = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if k: d[k] = ''.join(buf)
    return d
cand = read_fasta(f'{R}/out/rf_candidates.faa'); prots = {}
SHORT = {'GCA_056209055.1': 'T', 'GCA_934647165.1': 'C1', 'GCA_030529995.1': 'C2', 'GCA_963656345.1': 'C3', 'GCA_030526225.1': 'J', 'GCA_045599485.1': 'S'}
for l in list(open(f'{R}/out/rf_hmm.tsv'))[1:]:
    f = l.rstrip('\n').split('\t')
    if f[1] == 'std' and f[12] in ('RF1', 'RF2') and int(f[6]) >= 300 and f[0] in SHORT:
        prots[f'{SHORT[f[0]]}_{f[12]}'] = cand[f'{f[0]}|std|{f[2]}|{f[4]}-{f[5]}|{f[3]}']
for k, v in read_fasta(f'{R}/data/rf_queries.faa').items():
    if 'prfA' in k: prots[k.split('|')[0].replace('_prfA', '_RF1')] = v
    if 'prfB' in k: prots[k.split('|')[0].replace('_prfB', '_RF2')] = v
fa = f'{R}/out/rf/rf_proteins.faa'
with open(fa, 'w') as o:
    for k, v in prots.items(): o.write(f'>{k}\n{v}\n')
out = subprocess.run(['nice', '-n', '10', 'diamond', 'blastp', '-q', fa, '-d', f'{R}/data/sprot', '--threads', '2', '--quiet', '-k', '50', '-e', '1e-5',
                      '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'bitscore', 'stitle'], capture_output=True, text=True, check=True).stdout
hits = {}
for l in out.splitlines():
    q, sid, pid, ln, ql, sl, bits, st = l.split('\t'); hits.setdefault(q, []).append((float(bits), float(pid), int(ln), sid.split('|')[2], st.split(' OS=')[0].split(' ', 1)[1]))
print('## Swiss-Prot reciprocal check (DIAMOND, top 3; best RF1 = "chain release factor 1", best RF2 = "chain release factor 2")')
for k in prots:
    h = hits.get(k, [])
    b1 = next((x for x in h if 'release factor 1' in x[4] and 'mitochondrial' not in x[4].lower()), None)
    b2 = next((x for x in h if 'release factor 2' in x[4] and 'mitochondrial' not in x[4].lower()), None)
    f1 = f'  || best RF1 {b1[3]} {b1[1]:.0f}% ({b1[0]:.0f} bits, aln {b1[2]})' if b1 else '  || best RF1 none in top 50'
    f2 = f'  best RF2 {b2[3]} {b2[1]:.0f}% ({b2[0]:.0f} bits)' if b2 else '  best RF2 none in top 50'
    print(f'{k:8s} len={len(prots[k])} top: ' + ' | '.join(f'{x[3]} {x[1]:.0f}% {x[0]:.0f}b' for x in h[:3]) + f1 + f2)
al = Align.PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62'); al.open_gap_score = -11; al.extend_gap_score = -1; al.mode = 'local'
def pid(a, b):
    x = al.align(a, b)[0]; m = n = 0
    for (s1, e1), (s2, e2) in zip(*x.aligned):
        for i in range(e1 - s1): n += 1; m += a[s1 + i] == b[s2 + i]
    return m / max(1, n), n / min(len(a), len(b))
ks = list(prots)
print('\n## pairwise identity (identity over aligned columns / coverage of shorter)')
print('         ' + ' '.join(f'{k:>11s}' for k in ks))
for a in ks:
    print(f'{a:8s} ' + ' '.join(f'{pid(prots[a], prots[b])[0]:>6.0%}/{pid(prots[a], prots[b])[1]:>3.0%}' if a != b else f'{"-":>11s}' for b in ks))
def mapwin(ref, q, pos, w):
    x = al.align(ref, q)[0]; m = {}
    for (s1, e1), (s2, e2) in zip(*x.aligned):
        for i in range(e1 - s1): m[s1 + i + 1] = q[s2 + i]
    # insertions in q relative to ref are not shown; deletions shown as '-'
    return ''.join(m.get(p, '-') for p in range(pos - w, pos + w + 3))
ec1, ec2 = prots['Ec_RF1'], prots['Ec_RF2']
i_pxt = ec1.find('PAT', 150) + 1; i_spf = ec2.find('SPF') + 1; i_ggq = ec1.find('GGQ') + 1; i_ggq2 = ec2.find('GGQ') + 1
print(f'\n## motif windows (aligned to E. coli RF1 PAT at {i_pxt} +-6, E. coli RF2 SPF at {i_spf} +-6, GGQ +-4); "-" = deleted vs E. coli')
for k in ks:
    print(f'{k:8s} vsEcRF1 PxT-win {mapwin(ec1, prots[k], i_pxt, 6)}  vsEcRF2 SPF-win {mapwin(ec2, prots[k], i_spf, 6)}  GGQ-win {mapwin(ec1, prots[k], i_ggq, 4)}')
# raw context around the PxT region in each protein (substring search on the aligned window anchor)
print('\n## local alignment of T_RF1 with C2_RF1 and C3_RF1 (full)')
for o in ('C2_RF1', 'C3_RF1', 'Ec_RF1'):
    if o in prots and 'T_RF1' in prots: print(o); print(al.align(prots['T_RF1'], prots[o])[0])
