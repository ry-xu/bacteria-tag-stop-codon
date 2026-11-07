"""Rebuild the recheck/recheck2 conserved-protein query set (E. coli K-12 + B. subtilis 168, ~457 proteins
selected by gene name) without downloading: query IDs (org|gene|protein_id) are taken from recheck2's
tblastn outputs, full-length sequences from recheck2/out/rf/ecbs_prot.faa (complete RefSeq proteomes).
Writes data/queries.faa (conserved set) and data/rf_queries.faa (Ec/Bs prfA, prfB, prmC, Ec arfB).
Usage: python3 -I extract_queries.py <recheck3_dir> <recheck2_dir>"""
import sys
R, R2 = sys.argv[1], sys.argv[2]
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
prot = {k.split('|')[1]: v for k, v in read_fasta(f'{R2}/out/rf/ecbs_prot.faa').items()}
ids = []
for acc in ('GCF_000005845.2', 'GCF_000009045.1'):
    for line in open(f'{R2}/out/{acc}.tblastn.tsv'):
        q = line.split('\t', 1)[0]
        if q not in ids: ids.append(q)
nq = nr = 0
with open(f'{R}/data/queries.faa', 'w') as o, open(f'{R}/data/rf_queries.faa', 'w') as orf:
    for q in ids:
        pid = q.split('|')[-1]
        if q.startswith(('Ec|', 'Bs|')): o.write(f'>{q}\n{prot[pid]}\n'); nq += 1
        else: orf.write(f'>{q}\n{prot[pid]}\n'); nr += 1
print(nq, 'conserved queries;', nr, 'RF queries')
