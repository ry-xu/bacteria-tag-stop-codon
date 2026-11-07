"""prfA neighbourhood. pyrodigal (single mode, trained per genome, table 11) on each control genome that has an HMM RF1 call;
take K genes each side of prfA; annotate by DIAMOND blastp vs Swiss-Prot (top hit); tblastn each neighbour vs every target genome
(e<1e-10) and report the best hit (contig:start-end strand, %id, qcov). Writes out/rf/prfA_neigh_genes.faa and
out/rf/prfA_synteny.tsv. Then, per target, lists the contigs that carry the neighbours with the neighbour order along
the contig, the contig length and whether the run of neighbours stops at a contig end.
Usage: uv run --no-project --with pyrodigal python -I prfA_synteny.py <dir> K"""
import sys, subprocess, collections, os
import pyrodigal
R = sys.argv[1]; K = int(sys.argv[2])
targets = ['GCA_937907055.1', 'GCA_045628485.1', 'GCA_030526225.1', 'GCA_934647165.1', 'GCA_017577285.1', 'GCA_937886445.1', 'GCA_045617125.1', 'GCA_030529995.1']
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
loci = {}
for l in list(open(f'{R}/out/rf_hmm.tsv'))[1:]:
    f = l.rstrip('\n').split('\t')
    if f[1] == 'std' and f[12] == 'RF1' and f[0].startswith('GCA_'): loci[f[0]] = (f[2], int(f[4]), int(f[5]))
os.makedirs(f'{R}/out/rf', exist_ok=True)
nb = f'{R}/out/rf/prfA_neigh_genes.faa'; order = []
with open(nb, 'w') as o:
    for acc, (c, s, e) in loci.items():
        g = read_fasta(f'{R}/data/{acc}.fna')
        gf = pyrodigal.GeneFinder(); gf.train(*[x.encode() for x in g.values() if len(x) > 20000][:50])
        genes = gf.find_genes(g[c].encode())
        gl = [(x.begin, x.end, x.strand, x.translate().rstrip('*')) for x in genes]
        idx = max(range(len(gl)), key=lambda i: min(gl[i][1], e) - max(gl[i][0], s))
        for j in range(max(0, idx - K), min(len(gl), idx + K + 1)):
            b, en, st, p = gl[j]
            name = f'{acc}|{c}|{b}-{en}|{st:+d}|{j-idx:+d}'
            o.write(f'>{name}\n{p}\n'); order.append(name)
        print(acc, c, 'len', len(g[c]), 'genes on contig', len(gl), 'prfA index', idx, file=sys.stderr)
ann = {}
out = subprocess.run(['nice', '-n', '10', 'diamond', 'blastp', '-q', nb, '-d', f'{R}/data/sprot', '--threads', '3', '--quiet', '-k', '1', '-e', '1e-5',
                      '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'evalue', 'stitle'], capture_output=True, text=True, check=True).stdout
for l in out.splitlines():
    q, sid, pid, ev, st = l.split('\t')
    if q not in ann: ann[q] = f'{st.split(" OS=")[0].split(" ", 1)[1][:45]} [{sid.split("|")[2]} {float(pid):.0f}%]'
hits = {}
for t in targets:
    out = subprocess.run(['nice', '-n', '10', 'tblastn', '-query', nb, '-db', f'{R}/data/{t}', '-evalue', '1e-10', '-max_target_seqs', '5', '-num_threads', '3',
                          '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'], capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        f = l.split('\t'); q = f[0]
        if (q, t) in hits and hits[(q, t)][5] >= float(f[9]): continue
        sid = f[1].split('|')[1] if '|' in f[1] else f[1]
        hits[(q, t)] = (sid, int(f[6]), int(f[7]), float(f[2]), (int(f[5]) - int(f[4]) + 1) / int(f[10]), float(f[9]), int(f[11]))
with open(f'{R}/out/rf/prfA_synteny.tsv', 'w') as o:
    o.write('control_gene\tannotation\t' + '\t'.join(targets) + '\n')
    for q in order:
        cells = []
        for t in targets:
            h = hits.get((q, t)); cells.append('-' if h is None else f'{h[0]}:{h[1]}-{h[2]}({h[6]}) {h[3]:.0f}% c{h[4]:.0%}')
        o.write(f'{q}\t{ann.get(q, "no Swiss-Prot hit")}\t' + '\t'.join(cells) + '\n')
print('done', file=sys.stderr)
