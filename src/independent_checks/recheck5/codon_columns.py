"""Codon-to-residue alignment columns between a test genome T and another genome O (no gene calls).
T: six-frame ORFs with only TAA/TGA as stops (>=100 codons, stop-to-stop), translated with TAG -> X
(X is neutral in scoring, so the residue aligned to a TAG is not pulled toward any amino acid).
tblastn (std code) of all T ORFs against O (e<1e-10). For every aligned column record the T codon, the O codon
(read from O's DNA) and O's residue. Filters:
  * frame competition: a column is kept only if its HSP has the highest bitscore among all HSPs covering
    each of its 3 T nucleotides (removes antisense/shadow-frame ORF alignments);
  * per T codon, the column from the best-scoring HSP is kept; win10_ident = identity of the +-10 flanking
    columns (centre excluded), cols_to_hsp_edge = distance in columns to the nearer HSP end;
  * column must sit inside an ungapped block (no gap in either sequence within +-3 columns), so alignment
    placement is not driven by the residue at that column.
Each column is also flagged if the T codon lies inside a 'census' HSP (conserved Ec/Bs protein, e<1e-10, qcov>=0.5).
Usage: python3 -I codon_columns.py <dir> <T_acc> <O_acc> [threads]
Writes <dir>/out/cols/<T>__<O>.tsv. (recheck3: added column O_stops_pm20 = stops in O's translated HSP within +-20
columns of the site, excluding the site; used to drop alignments to control non-coding/shadow frames.)"""
import sys, os, subprocess
R, T, O = sys.argv[1:4]; TH = sys.argv[4] if len(sys.argv) > 4 else '3'
CODE = {}; b = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for x in b:
    for y in b:
        for z in b: CODE[x+y+z] = aas[i]; i += 1
COMP = str.maketrans('ACGT', 'TGCA')
def rc(s): return s.translate(COMP)[::-1]
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
gT = read_fasta(f'{R}/data/{T}.fna'); gO = read_fasta(f'{R}/data/{O}.fna')
os.makedirs(f'{R}/out/cols', exist_ok=True)
orf_fa = f'{R}/out/cols/{T}.tagsense_orfs.faa'
orfs = {}
if True:
    with open(orf_fa, 'w') as o:
        n = 0
        for c, s in gT.items():
            L = len(s)
            for strand, seq in ((1, s), (-1, rc(s))):
                for fr in range(3):
                    prot = []; start = fr
                    for p in range(fr, L - 2, 3):
                        cod = seq[p:p+3]
                        aa = 'X' if cod == 'TAG' else CODE.get(cod, 'X')
                        if aa == '*':
                            if len(prot) >= 100:
                                name = f'o{n}'; n += 1; orfs[name] = (c, strand, start); o.write(f'>{name}\n{"".join(prot)}\n')
                            prot = []; start = p + 3
                        else: prot.append(aa)
                    if len(prot) >= 100:
                        name = f'o{n}'; n += 1; orfs[name] = (c, strand, start); o.write(f'>{name}\n{"".join(prot)}\n')
def tcodon(name, qpos):
    """genomic 0-based start (on + strand) and codon (oriented) of query position qpos (1-based)"""
    c, strand, start = orfs[name]; L = len(gT[c]); o = start + 3 * (qpos - 1)
    if strand == 1: return c, strand, o, gT[c][o:o+3]
    g0 = L - o - 3; return c, strand, g0, rc(gT[c][g0:g0+3])
# census HSP spans in T (from the conserved-protein tblastn)
census = {}
for line in open(f'{R}/out/{T}.tblastn.tsv'):
    f = line.rstrip('\n').split('\t'); q = f[0]
    if not (q.startswith('Ec|') or q.startswith('Bs|')): continue
    qs, qe, ss, se, ev, qlen = int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[8]), int(f[10])
    if ev >= 1e-10 or (qe - qs + 1) / qlen < 0.5: continue
    s = f[1].split('|')[1] if '|' in f[1] else f[1]
    census.setdefault(s, []).append((min(ss, se) - 1, max(ss, se), 1 if se > ss else -1))
def in_census(c, g0, strand):
    # same strand AND same codon phase as a census HSP (frame-aware)
    for a, b, st in census.get(c, ()):
        if st != strand or not (a <= g0 and g0 + 3 <= b): continue
        if (st == 1 and (g0 - a) % 3 == 0) or (st == -1 and (b - (g0 + 3)) % 3 == 0): return True
    return False
bl = f'{R}/out/cols/{T}__{O}.tblastn.tsv'
if not os.path.exists(bl):
    subprocess.run(['nice', '-n', '10', 'tblastn', '-query', orf_fa, '-db', f'{R}/data/{O}', '-db_gencode', '1', '-seg', 'no', '-evalue', '1e-10',
                    '-max_target_seqs', '5', '-num_threads', TH, '-outfmt',
                    '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sframe qseq sseq', '-out', bl], check=True)
hsps = []
for line in open(bl):
    f = line.rstrip('\n').split('\t')
    hsps.append((f[0], f[1].split('|')[1] if '|' in f[1] else f[1], float(f[2]), int(f[4]), int(f[5]), int(f[6]), int(f[7]),
                 float(f[9]), int(f[10]), f[11], f[12]))
# frame competition: best HSP bitscore per T nucleotide
best = {c: {} for c in gT}
for h in hsps:
    name, _, _, qs, qe = h[:5]; bits = h[7]
    c, _, a, _ = tcodon(name, qs); _, _, z, _ = tcodon(name, qe)
    lo, hi = min(a, z), max(a, z) + 3
    bc = best[c]
    for p in range(lo, hi):
        if bc.get(p, 0) < bits: bc[p] = bits
cols = {}
for h in hsps:
    name, oc, pid, qs, qe, ss, se, bits, sframe, qseq, sseq = h
    oseq = gO[oc]; ostrand = 1 if sframe > 0 else -1
    qpos = qs; opos = ss
    n = len(qseq)
    gap = [qseq[k] == '-' or sseq[k] == '-' for k in range(n)]
    for k in range(n):
        qa, sa = qseq[k], sseq[k]
        if sa != '-':
            if ostrand == 1: ocod = oseq[opos-1:opos+2]; opos_cur = opos; opos += 3
            else: ocod = rc(oseq[opos-3:opos]); opos_cur = opos; opos -= 3
        if qa != '-' and sa != '-':
            if k >= 3 and k + 3 < n and not any(gap[k-3:k+4]):
                c, strand, g0, tcod = tcodon(name, qpos)
                bc = best[c]
                if all(bc.get(p, 0) <= bits for p in range(g0, g0 + 3)):
                    key = (c, g0, strand)
                    if key not in cols or cols[key][4] < bits:
                        assert (CODE.get(ocod, 'X') == sa) or ('N' in ocod) or sa == 'X', (ocod, sa)
                        assert tcod == 'TAG' and qa == 'X' or CODE.get(tcod, 'X') == qa or 'N' in tcod, (tcod, qa)
                        win = [j for j in range(max(0, k - 10), min(n, k + 11)) if j != k and qseq[j] != '-' and sseq[j] != '-']
                        wid = sum(1 for j in win if qseq[j] == sseq[j]) / max(1, len(win))
                        ostop20 = sseq[max(0, k - 20):k].count('*') + sseq[k + 1:k + 21].count('*')
                        cols[key] = (tcod, oc, opos_cur, ocod, bits, pid, sa, round(wid, 2), min(k, n - 1 - k), ostop20)
        if qa != '-': qpos += 1
with open(f'{R}/out/cols/{T}__{O}.tsv', 'w') as o:
    o.write('T_contig\tT_g0\tT_strand\tT_codon\tO_contig\tO_pos\tO_codon\tO_res\tbits\tpident\tin_census\twin10_ident\tcols_to_hsp_edge\tO_stops_pm20\n')
    for (c, g0, strand), v in sorted(cols.items()):
        o.write(f'{c}\t{g0}\t{strand}\t{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[6]}\t{v[4]}\t{v[5]}\t{int(in_census(c, g0, strand))}\t{v[7]}\t{v[8]}\t{v[9]}\n')
print(T, O, 'orfs', len(orfs), 'hsps', len(hsps), 'columns', len(cols), 'TAG columns', sum(1 for v in cols.values() if v[0] == 'TAG'))
