"""prfA neighbourhood: pyrodigal gene calls (table 11, single mode trained on the whole genome) in each
control genome; take the 8 genes on each side of prfA (the gene overlapping the HMM prfA call); annotate each by
blastp vs E. coli K-12 + B. subtilis 168 RefSeq proteomes; tblastn each against the two test genomes and
report best hit (contig, coords, strand, %id, qcov). Also prints, for each test genome, the gene order of hits
on the contig(s) carrying the neighbours.
Usage: uv run --with pyrodigal python -I prfA_synteny.py <recheck_dir>"""
import sys, subprocess, collections
import pyrodigal
R = sys.argv[1]; K = 8
loci = {'GCA_934796175.1': ('CAKVOV010000023.1', 30989, 32167),
        'GCA_004556155.1': ('SFIM01000024.1', 3755, 4933),
        'GCA_963624325.1': ('CAVFPH010000002.1', 73650, 74726)}
tests = ['GCA_024404275.1', 'GCA_045599485.1']
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
# reference proteome db for annotation
ref = f'{R}/out/rf/ecbs_prot'
with open(ref + '.faa', 'w') as o:
    for acc, org in (('GCF_000005845.2', 'Ec'), ('GCF_000009045.1', 'Bs')):
        for line in open(f'{R}/dl/{acc}/ncbi_dataset/data/{acc}/protein.faa'):
            if line.startswith('>'):
                pid, desc = line[1:].split(' ', 1); desc = desc.split(' [')[0].replace(' ', '_')
                o.write(f'>{org}|{pid}|{desc[:60]}\n')
            else: o.write(line)
subprocess.run(['makeblastdb', '-in', ref + '.faa', '-dbtype', 'prot', '-out', ref], check=True, capture_output=True)
nb = f'{R}/out/rf/prfA_neigh_genes.faa'
order = []
with open(nb, 'w') as o:
    for acc, (c, s, e) in loci.items():
        g = read_fasta(f'{R}/data/{acc}.fna')
        gf = pyrodigal.GeneFinder(); gf.train(b'TTAATTAATTAA'.join(x.encode() for x in g.values()))
        genes = gf.find_genes(g[c].encode())
        gl = [(x.begin, x.end, x.strand, x.translate().rstrip('*')) for x in genes]
        idx = max(range(len(gl)), key=lambda i: min(gl[i][1], e) - max(gl[i][0], s))
        for j in range(max(0, idx - K), min(len(gl), idx + K + 1)):
            b, en, st, p = gl[j]
            name = f'{acc}|{c}|{b}-{en}|{st:+d}|{j-idx:+d}' + ('|prfA' if j == idx else '')
            o.write(f'>{name}\n{p}\n'); order.append(name)
ann = {}
out = subprocess.run(['blastp', '-query', nb, '-db', ref, '-evalue', '1e-5', '-max_target_seqs', '1', '-outfmt', '6 qseqid sseqid pident evalue',
                      '-num_threads', '2'], capture_output=True, text=True, check=True).stdout
for l in out.splitlines():
    q, sid, pid, ev = l.split('\t')
    if q not in ann: ann[q] = f'{sid.split("|")[2][:40]}({float(pid):.0f}%)'
hits = {}
for t in tests:
    out = subprocess.run(['tblastn', '-query', nb, '-db', f'{R}/data/{t}', '-evalue', '1e-10', '-max_target_seqs', '3', '-num_threads', '2',
                          '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen'], capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        f = l.split('\t'); q = f[0]
        if (q, t) in hits and hits[(q, t)][-1] >= float(f[9]): continue
        sid = f[1].split('|')[1] if '|' in f[1] else f[1]
        qcov = (int(f[5]) - int(f[4]) + 1) / int(f[10])
        hits[(q, t)] = (sid, int(f[6]), int(f[7]), float(f[2]), qcov, float(f[9]))
print('control_gene\tannotation(best Ec/Bs hit)\t' + '\t'.join(tests))
for q in order:
    cells = []
    for t in tests:
        h = hits.get((q, t)); cells.append('-' if h is None else f'{h[0]}:{h[1]}-{h[2]} id={h[3]:.0f}% cov={h[4]:.0%}')
    print(f'{q}\t{ann.get(q, "no Ec/Bs hit")}\t' + '\t'.join(cells))
