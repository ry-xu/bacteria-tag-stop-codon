"""Stop-codon read-through test, same treatment for TAA/TAG/TGA. Genes: pyrodigal (std code) calls in out/genes/<acc>.genes.tsv.
For each complete gene of the test genome (stop codon present, >= 100 aa, stop >= 120 nt from the contig end) build
'gene protein + translation of the next 40 in-frame codons' (stops -> X). DIAMOND blastp (more-sensitive, e<1e-10) against the
pyrodigal proteomes of the reference genomes. Read-through call: the best hit covers the stop (query alignment starts >= 30 aa
before it) and extends >= 15 aa past it, the subject continues >= 15 aa past the aligned stop position, and >= 30% identity
over the downstream segment. Reported per stop codon of the test gene: genes with a hit covering the last 30 aa, read-through
count and fraction. Usage: uv run --no-project --with pyrodigal python -I readthrough.py <dir> test_acc ref1,ref2,..."""
import sys, subprocess, collections, os
R, T, refs = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
COMP = str.maketrans('ACGT', 'TGCA'); rc = lambda s: s.translate(COMP)[::-1]
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
def genes(acc):
    for l in open(f'{R}/out/genes/{acc}.genes.tsv'):
        c, b_, e_, st = l.rstrip('\n').split('\t'); yield c, int(b_), int(e_), int(st)
def tr(s): return ''.join('X' if CODE.get(s[i:i+3], 'X') == '*' else CODE.get(s[i:i+3], 'X') for i in range(0, len(s) - 2, 3))
os.makedirs(f'{R}/out/rt', exist_ok=True)
# reference proteomes
db = f'{R}/out/rt/refs_{"_".join(r[4:13] for r in refs)}'
with open(db + '.faa', 'w') as o:
    for r in refs:
        g = read_fasta(f'{R}/data/{r}.fna')
        for n, (c, b_, e_, st) in enumerate(genes(r)):
            s = g[c][b_-1:e_]; s = s if st == 1 else rc(s); p = ''.join(CODE.get(s[i:i+3], 'X') for i in range(0, len(s) - 2, 3)).rstrip('*')
            o.write(f'>{r}|{n}\n{p.replace("*", "X")}\n')
subprocess.run(['diamond', 'makedb', '--in', db + '.faa', '-d', db, '--quiet'], check=True)
g = read_fasta(f'{R}/data/{T}.fna'); q = f'{R}/out/rt/{T}.ext.faa'; meta = {}
with open(q, 'w') as o:
    for n, (c, b_, e_, st) in enumerate(genes(T)):
        L = len(g[c])
        if st == 1:
            if e_ + 120 > L: continue
            s = g[c][b_-1:e_ + 120]
        else:
            if b_ - 121 < 0: continue
            s = rc(g[c][b_ - 121:e_])
        gl = (e_ - b_ + 1) // 3; stop = s[3 * gl - 3:3 * gl]
        if stop not in ('TAA', 'TAG', 'TGA') or gl < 101: continue
        p = tr(s); meta[f'g{n}'] = (c, b_, e_, st, stop, gl - 1)   # gl-1 = protein length without stop; stop at position gl (1-based)
        o.write(f'>g{n}\n{p}\n')
out = subprocess.run(['nice', '-n', '10', 'diamond', 'blastp', '-q', q, '-d', db, '--more-sensitive', '--threads', '3', '--quiet', '-k', '1', '-e', '1e-10',
                      '--outfmt', '6', 'qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'slen', 'qseq_gapped', 'sseq_gapped'],   # audit fix: qseq/sseq are ungapped in DIAMOND
                     capture_output=True, text=True, check=True).stdout
res = collections.defaultdict(lambda: [0, 0]); rt_list = []
for l in out.splitlines():
    qid, sid, qs, qe, ss, se, sl, qseq, sseq = l.split('\t'); qs, qe, ss, se, sl = map(int, (qs, qe, ss, se, sl))
    c, b_, e_, st, stop, plen = meta[qid]; sp = plen + 1      # stop position in query
    if qe < plen - 30 or qs > plen - 30: continue            # hit must cover the C-terminal 30 aa
    res[stop][0] += 1
    # walk alignment to find the subject position aligned to the stop and identity downstream
    qp, sp_ = qs, ss; down = ident = 0; s_at_stop = None
    for x, y in zip(qseq, sseq):
        if x != '-' and qp == sp: s_at_stop = sp_
        if x != '-' and y != '-' and qp > sp: down += 1; ident += x == y
        if x != '-': qp += 1
        if y != '-': sp_ += 1
    if qe >= sp + 15 and s_at_stop is not None and sl - s_at_stop >= 15 and down >= 15 and ident / down >= 0.3:
        res[stop][1] += 1; rt_list.append((stop, c, b_, e_, st, plen, qe - sp, round(ident / down, 2), sid))
print(f'## {T} vs {",".join(refs)}: complete genes with a reference hit covering their last 30 aa, and read-through calls')
for s in ('TAA', 'TAG', 'TGA'):
    n, k = res[s]; print(f'   {s}: genes={n} read-through={k} ({k/max(1,n):.1%})')
with open(f'{R}/out/rt/{T}.readthrough.tsv', 'w') as o:
    o.write('stop\tcontig\tbegin\tend\tstrand\tprot_len\taa_aligned_past_stop\tdownstream_identity\tref_hit\n')
    for r in rt_list: o.write('\t'.join(map(str, r)) + '\n')
