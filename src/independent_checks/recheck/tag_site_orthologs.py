"""For each in-alignment TAG site in a test genome, translate a +-50-codon window (TAG->X; window
trimmed at TAA/TGA), tblastn it against relative genomes, and report the residue aligned to the X.
Usage: python -I tag_site_orthologs.py <recheck_dir> <test_acc> <rel_acc> [<rel_acc> ...]"""
import sys, subprocess, collections, os
R, test = sys.argv[1], sys.argv[2]; rels = sys.argv[3:]
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
CODE['TAG'] = 'X'
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
g = read_fasta(f'{R}/data/{test}.fna')
sites = [l.rstrip('\n').split('\t') for l in open(f'{R}/out/{test}.stop_sites.tsv')][1:]
sites = [s for s in sites if s[3] == 'TAG']
qf = f'{R}/out/{test}.tagwin.faa'; centers = {}
with open(qf, 'w') as o:
    for s in sites:
        c, strand, p = s[0], int(s[1]), int(s[2]); seq = g[c]
        if strand == 1: cods = [seq[p-1+3*k:p+2+3*k] for k in range(-50, 51) if p-1+3*k >= 0 and p+2+3*k <= len(seq)]; off = min(50, (p-1)//3)
        else: cods = [seq[p-3-3*k:p-3*k].translate(COMP)[::-1] for k in range(-50, 51) if p-3-3*k >= 0 and p-3*k <= len(seq)]; off = min(50, (p-3)//3)
        prot = ''.join(CODE.get(x, 'X') if len(x) == 3 else '' for x in cods)
        assert prot[off] == 'X', (s, prot, off)
        lo = prot.rfind('*', 0, off) + 1; hi = prot.find('*', off); hi = len(prot) if hi < 0 else hi
        name = f'{c}:{strand}:{p}'; centers[name] = off - lo + 1  # 1-based
        o.write(f'>{name}\n{prot[lo:hi]}\n')
res = {}
for rel in rels:
    out = subprocess.run(['tblastn', '-query', qf, '-db', f'{R}/data/{rel}', '-evalue', '1e-5', '-max_target_seqs', '5',
                          '-num_threads', '3', '-outfmt', '6 qseqid sseqid qstart qend evalue bitscore qseq sseq'],
                         capture_output=True, text=True, check=True).stdout
    best = {}
    for line in out.splitlines():
        q, ss, qs, qe, ev, bits, qseq, sseq = line.split('\t'); qs, qe, bits = int(qs), int(qe), float(bits)
        cpos = centers[q]
        if not (qs <= cpos <= qe): continue
        if q in best and best[q][0] >= bits: continue
        k = qs - 1
        for qa, sa in zip(qseq, sseq):
            if qa != '-':
                k += 1
                if k == cpos: best[q] = (bits, sa); break
    for q in centers: res[(q, rel)] = best.get(q, (0, '.'))[1]
with open(f'{R}/out/{test}.tag_ortholog_residues.tsv', 'w') as o:
    o.write('site\tbest_query\tquery_pos\tquery_residue\t' + '\t'.join(rels) + '\n')
    for s in sites:
        name = f'{s[0]}:{s[1]}:{s[2]}'
        o.write(f'{name}\t{s[4]}\t{s[5]}\t{s[6]}\t' + '\t'.join(res[(name, r)] for r in rels) + '\n')
print(test)
for r in rels:
    cnt = collections.Counter(res[(q, r)] for q in centers)
    print(' ', r, 'n_aligned=', sum(v for k, v in cnt.items() if k != '.'), dict(cnt.most_common()))
