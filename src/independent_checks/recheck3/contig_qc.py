"""Per-contig QC of test genomes: length, GC, SPAdes coverage (if in header), fraction of the contig aligned (megablast,
e<1e-20) to a reference genome set and mean identity, and in-alignment TAG count (from step 1). Summarises TAG-bearing
contigs vs the rest, and lists every TAG-bearing contig.
Usage: python3 -I contig_qc.py <dir> test_acc ref_acc[,ref_acc...]"""
import sys, subprocess, re, collections, statistics
R, T, refs = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
def read_fasta(p):
    d = {}; h = {}; k = None; buf = []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if k: d[k] = ''.join(buf).upper()
            k = line[1:].split()[0]; h[k] = line; buf = []
        else: buf.append(line)
    if k: d[k] = ''.join(buf).upper()
    return d, h
g, hdr = read_fasta(f'{R}/data/{T}.fna')
tag = collections.Counter(l.split('\t')[0] for l in list(open(f'{R}/out/{T}.stop_sites.tsv'))[1:] if l.split('\t')[3] == 'TAG')
cov = {}
for c, h in hdr.items():
    m = re.search(r'_cov_([0-9.]+)', h); cov[c] = float(m.group(1)) if m else None
res = {}
for ref in refs:
    out = subprocess.run(['nice', '-n', '10', 'blastn', '-task', 'megablast', '-query', f'{R}/data/{T}.fna', '-db', f'{R}/data/{ref}', '-evalue', '1e-20',
                          '-num_threads', '1', '-outfmt', '6 qseqid pident length qstart qend', '-max_target_seqs', '5', '-max_hsps', '50'],
                         capture_output=True, text=True, check=True).stdout
    cv = collections.defaultdict(lambda: None)
    per = collections.defaultdict(list)
    for l in out.splitlines():
        q, pid, ln, a, b = l.split('\t'); per[q].append((int(a), int(b), float(pid)))
    for c in g:
        mask = bytearray(len(g[c])); ids = []
        for a, b, pid in sorted(per[c], key=lambda x: -x[2]):
            new = sum(1 for i in range(a - 1, b) if not mask[i])
            if new > 0.5 * (b - a + 1):
                for i in range(a - 1, b): mask[i] = 1
                ids.append((pid, b - a + 1))
        af = sum(mask) / len(g[c]); mid = sum(p * n for p, n in ids) / max(1, sum(n for _, n in ids))
        res.setdefault(c, {})[ref] = (af, mid)
def row(c):
    gc = (g[c].count('G') + g[c].count('C')) / len(g[c])
    best = max(res[c].items(), key=lambda x: x[1][0])
    return gc, best
print(f'## {T} vs {refs}')
for lab, cs in (('TAG-bearing', [c for c in g if tag[c]]), ('other', [c for c in g if not tag[c]])):
    L = sum(len(g[c]) for c in cs)
    af = sum(max(v[0] for v in res[c].values()) * len(g[c]) for c in cs) / L
    cvs = [cov[c] for c in cs if cov[c] is not None]
    print(f'{lab}: contigs={len(cs)} bp={L} GC={sum(row(c)[0]*len(g[c]) for c in cs)/L:.3f} aligned_frac_to_best_ref={af:.2f} '
          + (f'cov median={statistics.median(cvs):.1f}' if cvs else ''))
allcov = [cov[c] for c in g if cov[c] is not None]
if allcov: print(f'all-contig coverage quartiles: {statistics.quantiles(allcov, n=4)}')
print('contig\tlen\tGC\tcov\tTAG\t' + '\t'.join(f'{r}:aligned_frac/id' for r in refs))
for c in sorted(g, key=lambda c: -tag[c]):
    if not tag[c]: continue
    gc = row(c)[0]
    print(f'{c}\t{len(g[c])}\t{gc:.3f}\t{cov[c]}\t{tag[c]}\t' + '\t'.join(f'{res[c][r][0]:.2f}/{res[c][r][1]:.1f}' for r in refs))
