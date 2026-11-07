"""Search raw DNA for RF1 fragments/pseudogenes. Queries: full-length RF1 and RF2 ORFs called by rf_hmm.py in the
controls (std mode, call RF1/RF2) + E. coli/B. subtilis prfA/prfB. tblastn (seg off, e<=10, std code and also
code with TAG read as sense is unnecessary: tblastn shows TAG as '*' inside HSPs anyway) vs each target genome.
Hits are merged into loci (same contig, overlapping +-300 nt); per locus the best RF1-query and best RF2-query
bitscores are reported -> class of the locus. Any locus with RF1 bits > RF2 bits, or any RF-like locus not
overlapping a full-length RF2 call, is listed.
Usage: python3 -I rf_fragments.py <dir> target_acc..."""
import sys, subprocess, collections
R = sys.argv[1]; targets = sys.argv[2:]
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
calls = {}
for l in list(open(f'{R}/out/rf_hmm.tsv'))[1:]:
    f = l.rstrip('\n').split('\t')
    if f[1] == 'std' and f[12] in ('RF1', 'RF2') and int(f[6]) >= 300: calls[f'{f[0]}|std|{f[2]}|{f[4]}-{f[5]}|{f[3]}'] = f[12]
cand = read_fasta(f'{R}/out/rf_candidates.faa')
q = f'{R}/out/rf/frag_queries.faa'
import os; os.makedirs(f'{R}/out/rf', exist_ok=True)
with open(q, 'w') as o:
    for k, cl in calls.items():
        if k.startswith('GCF_') or k.split('|')[0] in targets and cl == 'RF1': pass
        if k.split('|')[0] in targets: continue   # do not use a target's own proteins as queries
        o.write(f'>{cl}|{k.split("|")[0]}|{k.split("|")[2]}\n{cand[k]}\n')
    for k, v in read_fasta(f'{R}/data/rf_queries.faa').items():
        if 'prfA' in k: o.write(f'>RF1|{k}\n{v}\n')
        elif 'prfB' in k: o.write(f'>RF2|{k}\n{v}\n')
for t in targets:
    out = subprocess.run(['nice', '-n', '10', 'tblastn', '-query', q, '-db', f'{R}/data/{t}', '-seg', 'no', '-evalue', '10', '-max_target_seqs', '50',
                          '-num_threads', '1', '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen sseq'],
                         capture_output=True, text=True, check=True).stdout
    hits = []
    for l in out.splitlines():
        f = l.split('\t'); s = f[1].split('|')[1] if '|' in f[1] else f[1]
        a, b = sorted((int(f[6]), int(f[7])))
        hits.append((s, a, b, f[0], float(f[2]), int(f[4]), int(f[5]), float(f[8]), float(f[9]), f[11].count('*')))
    hits.sort(); loci = []
    for h in hits:
        if loci and loci[-1][0] == h[0] and h[1] <= loci[-1][2] + 300: loci[-1][2] = max(loci[-1][2], h[2]); loci[-1][3].append(h)
        else: loci.append([h[0], h[1], h[2], [h]])
    print(f'== {t}: {len(loci)} loci with any RF1/RF2 query hit (e<=10)')
    for s, a, b, hs in loci:
        b1 = max([h for h in hs if h[3].startswith('RF1')], key=lambda h: h[8], default=None)
        b2 = max([h for h in hs if h[3].startswith('RF2')], key=lambda h: h[8], default=None)
        if min(h[7] for h in hs) > 1e-2: continue   # recheck4 fix: was max() (skipped any locus that also had a weak hit)
        fmt = lambda h: '-' if h is None else f'{h[3].split("|")[1]} q{h[5]}-{h[6]} id={h[4]:.0f}% e={h[7]:.1e} bits={h[8]:.0f} stops={h[9]}'
        cls = 'RF1>RF2' if (b1 and (not b2 or b1[8] > b2[8])) else 'RF2>=RF1'
        print(f'  {s}:{a}-{b} {cls} | best RF1 query: {fmt(b1)} | best RF2 query: {fmt(b2)}')
