"""TAG columns of a test genome T that align to a stop codon ('*') in a control O (out/cols tables): these are TAGs at
the homologous position of the control's gene end. For each, distance (codons) in T from the TAG to the next in-frame
TAA/TGA and the number of further TAGs before it; compared with the same statistic for T's TAA/TGA columns aligned to
control stops is meaningless (they end the ORF), so compare to random in-frame positions: P(next TAA/TGA within k codons).
Usage: python3 -I tag_at_ctrl_stop.py <dir> T O[,O...]"""
import sys, collections
R, T, Os = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
COMP = str.maketrans('ACGT', 'TGCA')
def read_fasta(p):
    d = {}; k = None; buf = []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if k: d[k] = ''.join(buf).upper()
            k = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if k: d[k] = ''.join(buf).upper()
    return d
g = read_fasta(f'{R}/data/{T}.fna')
def downstream(c, g0, strand, n=200):
    s = g[c]
    if strand == '1': return [s[g0 + 3 * i: g0 + 3 * i + 3] for i in range(1, n) if g0 + 3 * i + 3 <= len(s)]
    out = []
    for i in range(1, n):
        a = g0 - 3 * i
        if a < 0: break
        out.append(s[a:a + 3].translate(COMP)[::-1])
    return out
sites = {}; sense = {}
for O in Os:
    for l in list(open(f'{R}/out/cols/{T}__{O}.tsv'))[1:]:
        r = l.rstrip('\n').split('\t')
        if r[3] != 'TAG': continue
        k = (r[0], int(r[1]), r[2])
        if r[7] == '*': sites.setdefault(k, set()).add(r[6])
        else: sense.setdefault(k, set()).add(r[7])
stop_only = {k: v for k, v in sites.items() if k not in sense}
dist = []
for (c, g0, st), ocods in sorted(stop_only.items()):
    d = downstream(c, g0, st)
    nx = next((i for i, x in enumerate(d) if x in ('TAA', 'TGA')), None)
    dist.append(nx)
    print(f'   {c}:{g0}:{st} control stop codon(s) {",".join(sorted(ocods))} -> next TAA/TGA in T after {"none<200" if nx is None else nx + 1} codons; TAGs before it: {d[:nx if nx is not None else 200].count("TAG")}')
dd = [x + 1 for x in dist if x is not None]
print(f'{T}: TAG columns aligned only to control stops: {len(stop_only)}; (also aligned to a sense residue in another control: {len(set(sites) - set(stop_only))})')
if dd:
    dd.sort(); print(f'   codons to next TAA/TGA: median {dd[len(dd)//2]}, <=3: {sum(1 for x in dd if x <= 3)}, <=10: {sum(1 for x in dd if x <= 10)}, n={len(dd)}')
# background: from random TAG-free in-frame positions inside T's ORFs? use all T sense TAG sites' downstream as reference
ref = []
for (c, g0, st) in sense:
    d = downstream(c, g0, st); nx = next((i for i, x in enumerate(d) if x in ('TAA', 'TGA')), None)
    if nx is not None: ref.append(nx + 1)
ref.sort()
if ref: print(f'   reference (TAG columns aligned to a sense residue): n={len(ref)} median {ref[len(ref)//2]}, <=3: {sum(1 for x in ref if x <= 3)}, <=10: {sum(1 for x in ref if x <= 10)}')
