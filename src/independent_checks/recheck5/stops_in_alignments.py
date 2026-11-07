"""(Copied unchanged in logic from recheck/scripts.) Count stop codons that tblastn aligned INSIDE high-scoring HSPs of conserved proteins.
Usage: python -I stops_in_alignments.py <recheck_dir> <acc> [<acc> ...]
Filters: evalue < 1e-10, query coverage (qend-qstart+1)/qlen >= 0.5, query is one of the conserved set (Ec|/Bs| ids).
Codons are deduplicated by genomic position (contig, strand, codon start), so overlapping
Ec/Bs hits or paralog hits to the same locus are counted once. Writes per-genome site tables."""
import sys, collections, os
R = sys.argv[1]; accs = sys.argv[2:]
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
STOPS = {'TAA', 'TAG', 'TGA'}

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

summary = []
for acc in accs:
    g = read_fasta(f'{R}/data/{acc}.fna')
    codons = {}      # key -> (codon, best bitscore, row info)
    hsps_used = 0; loci_q = collections.Counter()
    term = collections.Counter()   # first stop after C-terminal HSP end
    for line in open(f'{R}/out/{acc}.tblastn.tsv'):
        f = line.rstrip('\n').split('\t')
        q, s = f[0], (f[1].split('|')[1] if '|' in f[1] else f[1])
        if not (q.startswith('Ec|') or q.startswith('Bs|')): continue
        qs, qe, ss, se = int(f[4]), int(f[5]), int(f[6]), int(f[7])
        ev, bits, qlen, slen, sframe = float(f[8]), float(f[9]), int(f[10]), int(f[11]), int(f[12])
        qseq, sseq = f[13], f[14]
        if ev >= 1e-10 or (qe - qs + 1) / qlen < 0.5: continue
        hsps_used += 1
        seq = g[s]; strand = 1 if sframe > 0 else -1
        pos = ss; qpos = qs
        ntot = len(sseq) - sseq.count('-'); ci = 0
        for qa, sa in zip(qseq, sseq):
            if sa != '-':
                ci += 1
                if strand == 1:
                    cod = seq[pos-1:pos+2]; key = (s, 1, pos); pos += 3
                else:
                    cod = seq[pos-3:pos].translate(COMP)[::-1]; key = (s, -1, pos); pos -= 3
                if sa == '*': assert cod in STOPS, (acc, q, s, cod)
                elif cod in STOPS: raise SystemExit(f'non-* residue on stop codon {acc} {q} {s} {key} {cod} {sa}')
                prev = codons.get(key)
                if prev is None or bits > prev[1]:
                    codons[key] = (cod, bits, q, qpos if qa != '-' else None, qa, ev, qlen, qs, qe, min(ci - 1, ntot - ci))
            if qa != '-': qpos += 1
        # terminal stop: HSP ends within 10 aa of query C-terminus -> first in-frame stop downstream (<= 600 nt)
        if qe >= qlen - 10:
            p = pos  # next codon start after alignment end
            for _ in range(200):
                if strand == 1:
                    if p + 2 > len(seq): break
                    cod = seq[p-1:p+2]; p += 3
                else:
                    if p - 3 < 0: break
                    cod = seq[p-3:p].translate(COMP)[::-1]; p -= 3
                if cod in STOPS: term[(s, strand, p)] = cod; break
    # summarise
    n_cod = len(codons)
    st = collections.Counter(v[0] for v in codons.values() if v[0] in STOPS)
    termc = collections.Counter(term.values())
    # loci with >=1 in-alignment stop, by gene name of best hit
    genes_all = set(v[2].split('|')[1] for v in codons.values())
    genes_tag = set(v[2].split('|')[1] for v in codons.values() if v[0] == 'TAG')
    genes_taa = set(v[2].split('|')[1] for v in codons.values() if v[0] == 'TAA')
    genes_tga = set(v[2].split('|')[1] for v in codons.values() if v[0] == 'TGA')
    internal = collections.Counter(v[0] for v in codons.values() if v[0] in STOPS and v[9] >= 10)
    cu = collections.Counter(v[0] for v in codons.values())
    with open(f'{R}/out/{acc}.aligned_codon_usage.tsv', 'w') as o:
        o.write('codon\tcount\tper_1000\n')
        for c, n in cu.most_common(): o.write(f'{c}\t{n}\t{n/n_cod*1000:.3f}\n')
    row = dict(acc=acc, contigs=len(g), hsps=hsps_used, gene_families=len(genes_all), aligned_codons=n_cod,
               TAA=st['TAA'], TAG=st['TAG'], TGA=st['TGA'],
               TAA_per100k=round(st['TAA']/n_cod*1e5, 1), TAG_per100k=round(st['TAG']/n_cod*1e5, 1), TGA_per100k=round(st['TGA']/n_cod*1e5, 1),
               genes_with_TAA=len(genes_taa), genes_with_TAG=len(genes_tag), genes_with_TGA=len(genes_tga),
               internal10_TAA=internal['TAA'], internal10_TAG=internal['TAG'], internal10_TGA=internal['TGA'],
               term_TAA=termc['TAA'], term_TAG=termc['TAG'], term_TGA=termc['TGA'])
    summary.append(row)
    with open(f'{R}/out/{acc}.stop_sites.tsv', 'w') as o:
        o.write('contig\tstrand\tcodon_start\tcodon\tbest_query\tquery_pos\tquery_residue\tbitscore\tevalue\tqlen\tqstart\tqend\tcontig_len\tcodons_to_hsp_edge\n')
        for k, v in sorted(codons.items()):
            if v[0] in STOPS:
                o.write(f'{k[0]}\t{k[1]}\t{k[2]}\t{v[0]}\t{v[2]}\t{v[3]}\t{v[4]}\t{v[1]}\t{v[5]}\t{v[6]}\t{v[7]}\t{v[8]}\t{len(g[k[0]])}\t{v[9]}\n')
    # per-contig aligned codons
    pc = collections.Counter(k[0] for k in codons)
    with open(f'{R}/out/{acc}.contig_aligned.tsv', 'w') as o:
        o.write('contig\tcontig_len\taligned_codons\tTAA\tTAG\tTGA\n')
        sc = collections.defaultdict(collections.Counter)
        for k, v in codons.items():
            if v[0] in STOPS: sc[k[0]][v[0]] += 1
        for c, n in pc.most_common():
            o.write(f'{c}\t{len(g[c])}\t{n}\t{sc[c]["TAA"]}\t{sc[c]["TAG"]}\t{sc[c]["TGA"]}\n')
cols = list(summary[0].keys())
with open(f'{R}/out/summary_stops.tsv', 'w') as o:
    o.write('\t'.join(cols) + '\n')
    for r in summary: o.write('\t'.join(str(r[c]) for c in cols) + '\n')
print('\t'.join(cols))
for r in summary: print('\t'.join(str(r[c]) for c in cols))
