"""Redundant-segment map of a bin (two-strain / congeneric contamination signal): from blastn dc-megablast of the bin
against itself (out/dup/<acc>.self.tsv), keep non-self hits >= 500 nt with identity in [LO, HI); per contig report the
covered fraction and partner contigs. Prints totals and, for a list of contigs of interest, their redundancy.
Usage: python3 -I self_dup.py <dir> acc LO HI [contig ...]"""
import sys, collections
R, a, LO, HI = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]); want = sys.argv[5:]
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
g = read_fasta(f'{R}/data/{a}.fna')
cov = collections.defaultdict(set); part = collections.defaultdict(collections.Counter)
for l in open(f'{R}/out/dup/{a}.self.tsv'):
    f = l.split('\t'); q = f[0]; s = f[1].split('|')[1] if '|' in f[1] else f[1]
    pid = float(f[2]); L = int(f[3]); qa, qb = sorted((int(f[4]), int(f[5])))
    if q == s or L < 500 or not (LO <= pid < HI): continue
    cov[q].update(range(qa, qb + 1)); part[q][s] += L
tot = sum(len(v) for v in cov.values())
print(f'{a}: contigs with redundant segments ({LO}-{HI}% id, >=500 nt): {len(cov)}; covered nt {tot} ({tot/sum(map(len,g.values())):.1%} of bin)')
for c in sorted(cov, key=lambda c: -len(cov[c]))[:25]:
    print(f'   {c} len={len(g[c])} covered={len(cov[c])} ({len(cov[c])/len(g[c]):.0%}) partners ' + ' '.join(f'{s}:{n}' for s, n in part[c].most_common(3)))
if want:
    print('contigs of interest:')
    for c in want: print(f'   {c} len={len(g[c])} redundant={len(cov.get(c, ()))} ({len(cov.get(c, ()))/len(g[c]):.0%}) partners ' + ' '.join(f'{s}:{n}' for s, n in part[c].most_common(3)))
