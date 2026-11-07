"""Rough ANIb: split query genome into 1020-nt fragments, blastn vs subject genome db, keep best hit per
fragment with alignment >= 70% of fragment length and identity >= 30%; ANI = mean identity; AF = fraction aligned.
Usage: python3 -I ani.py <dir> query_acc subject_acc"""
import sys, subprocess
R, q, s = sys.argv[1:4]
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
g = read_fasta(f'{R}/data/{q}.fna'); fr = f'{R}/out/ani/{q}.frag.fna'
with open(fr, 'w') as o:
    for c, sq in g.items():
        for i in range(0, len(sq) - 1019, 1020): o.write(f'>{c}_{i}\n{sq[i:i+1020]}\n')
nf = sum(1 for l in open(fr) if l.startswith('>'))
out = subprocess.run(['nice', '-n', '10', 'blastn', '-task', 'blastn', '-query', fr, '-db', f'{R}/data/{s}', '-xdrop_gap', '150', '-penalty', '-1',
     '-reward', '1', '-evalue', '1e-15', '-max_target_seqs', '1', '-num_threads', '2', '-dust', 'no',
     '-outfmt', '6 qseqid pident length'], capture_output=True, text=True, check=True).stdout
best = {}
for l in out.splitlines():
    a, pid, ln = l.split('\t'); pid = float(pid); ln = int(ln)
    if ln >= 714 and pid >= 30 and a not in best: best[a] = pid
ids = list(best.values())
print(f'{q} vs {s}: ANI={sum(ids)/max(1,len(ids)):.2f}% fragments_aligned={len(ids)}/{nf} ({len(ids)/nf:.0%})')
