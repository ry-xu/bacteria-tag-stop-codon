"""Gene-call-free HMM search for class-I release factors.
Six-frame stop-to-stop ORFs (>=80 aa) in two codes: 'std' (TAA/TAG/TGA = stop) and
'tagsense' (only TAA/TGA = stop, TAG -> X). hmmsearch with PF03462 (PCRF), PF00472 (RF-1),
TIGR00019 (prfA), TIGR00020 (prfB). For each ORF with a PCRF or RF-1 hit, report all four scores,
and the residues aligned to E. coli RF1 P188-A189-T190 (PxT) and RF2 S205-P206-F207 (SPF), and GGQ.
Usage: uv run --with pyhmmer --with biopython python -I rf_hmm.py <recheck_dir> <hmm_dir> acc...
"""
import sys, collections
import pyhmmer
from pyhmmer.easel import TextSequence, Alphabet
from pyhmmer.plan7 import HMMFile
from Bio import Align
from Bio.Align import substitution_matrices
R, H = sys.argv[1], sys.argv[2]; accs = sys.argv[3:]
CODE = {}
bases = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODE[a+b+c] = aas[i]; i += 1
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
def orfs(g, mode, minlen=80):
    for c, s in g.items():
        L = len(s)
        for strand, seq in ((1, s), (-1, s.translate(COMP)[::-1])):
            for fr in range(3):
                prot = []; start = fr
                for p in range(fr, L - 2, 3):
                    cod = seq[p:p+3]
                    aa = CODE.get(cod, 'X')
                    if mode == 'tagsense' and cod == 'TAG': aa = 'X'
                    if aa == '*':
                        if len(prot) >= minlen: yield c, strand, start, p, ''.join(prot)
                        prot = []; start = p + 3
                    else: prot.append(aa)
                if len(prot) >= minlen: yield c, strand, start, (L - fr) // 3 * 3 + fr, ''.join(prot)
def coords(c, strand, start, end, L):
    # 1-based genomic coords of ORF (excluding stop)
    return (start + 1, end) if strand == 1 else (L - end + 1, L - start)
alpha = Alphabet.amino()
hmms = []
for f in ('PF03462.hmm', 'PF00472.hmm', 'TIGR00019.hmm', 'TIGR00020.hmm'):
    with HMMFile(f'{H}/{f}') as hf: hmms.extend(list(hf))
ref = {}
for line in open(f'{R}/data/rf_queries.faa'):
    line = line.strip()
    if line.startswith('>'): k = line[1:].split('|')[0]; ref[k] = ''
    else: ref[k] += line
aligner = Align.PairwiseAligner(); aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -11; aligner.extend_gap_score = -1; aligner.mode = 'local'
def map_motif(orf, refseq, positions):
    aln = aligner.align(refseq.replace('X', 'A'), orf)[0]
    m = {}
    for (rs, re_), (qs, qe) in zip(*aln.aligned):
        for k in range(re_ - rs): m[rs + k + 1] = orf[qs + k]
    return ''.join(m.get(p, '-') for p in positions), aln.score
out = open(f'{R}/out/rf_hmm.tsv', 'w')
out.write('genome\tmode\tcontig\tstrand\tstart\tend\taa_len\tPCRF\tRF-1\tTIGR00019_prfA\tTIGR00020_prfB\tcall\tEcRF1_188-190(PAT)\tEcRF2_205-207(SPF)\tGGQ(RF1_233-235)\tsw_vs_EcRF1\tsw_vs_EcRF2\tX_count\n')
fa = open(f'{R}/out/rf_candidates.faa', 'w')
for acc in accs:
    g = read_fasta(f'{R}/data/{acc}.fna')
    for mode in ('std', 'tagsense'):
        recs = []; seqs = []
        for i, (c, strand, st, en, prot) in enumerate(orfs(g, mode)):
            name = f'{acc}|{mode}|{i}'
            recs.append((c, strand, st, en, prot))
            seqs.append(TextSequence(name=name.encode(), sequence=prot).digitize(alpha))
        sb = pyhmmer.easel.DigitalSequenceBlock(alpha, seqs)
        score = collections.defaultdict(dict)
        for th in pyhmmer.hmmsearch(hmms, sb, cpus=3, E=1e-3):
            q = th.query
            hname = (q.name if hasattr(q, 'name') else th.query_name)
            hname = hname.decode() if isinstance(hname, bytes) else hname
            for h in th:
                n = h.name.decode() if isinstance(h.name, bytes) else h.name
                score[n][hname] = round(h.score, 1)
        for n, sc in score.items():
            if 'PCRF' not in sc and 'RF-1' not in sc: continue
            i = int(n.split('|')[2]); c, strand, st, en, prot = recs[i]
            s1, e1 = coords(c, strand, st, en, len(g[c]))
            a, b = sc.get('prfA', 0), sc.get('prfB', 0)
            call = 'RF1(prfA)' if a > b and a > 300 else 'RF2(prfB)' if b > a and b > 300 else 'other/partial'
            pxt, sw1 = map_motif(prot, ref['P0A7I0'], [188, 189, 190])
            spf, sw2 = map_motif(prot, ref['P07012'], [205, 206, 207])
            ggq, _ = map_motif(prot, ref['P0A7I0'], [233, 234, 235])
            out.write(f'{acc}\t{mode}\t{c}\t{strand}\t{s1}\t{e1}\t{len(prot)}\t{sc.get("PCRF","")}\t{sc.get("RF-1","")}\t{sc.get("prfA","")}\t{sc.get("prfB","")}\t{call}\t{pxt}\t{spf}\t{ggq}\t{sw1:.0f}\t{sw2:.0f}\t{prot.count("X")}\n')
            fa.write(f'>{acc}|{mode}|{c}|{s1}-{e1}|{strand}\n{prot}\n')
    print('done', acc, flush=True)
out.close(); fa.close()
