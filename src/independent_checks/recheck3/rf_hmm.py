"""Gene-call-free HMM search for class-I release factors RF1 (prfA) / RF2 (prfB).
HMMs: data/rf12_sprot.hmm (RF1_sprot, RF2_sprot; built locally from Swiss-Prot by build_rf_hmms.py; RF1-vs-RF2 score
margin over all 930 Swiss-Prot RF1/RF2 is >=132 bits in the right direction).
Six-frame stop-to-stop ORFs (>= MINLEN aa) in two codes: 'std' (TAA/TAG/TGA = stop) and 'tagsense' (TAG read as X).
For each ORF hit by either HMM (E<1e-3) report both scores, call (RF1 if RF1-RF2 >= 50 bits, RF2 if <= -50, else
ambiguous), the residues aligned to E. coli RF1 PxT (P-A-T) and RF2 SPF, GGQ, and SW scores vs E. coli RF1/RF2.
Usage: uv run --no-project --with pyhmmer --with biopython python -I rf_hmm.py <dir> MINLEN acc..."""
import sys, collections
import pyhmmer
from pyhmmer.easel import TextSequence, Alphabet, DigitalSequenceBlock
from pyhmmer.plan7 import HMMFile
from Bio import Align
from Bio.Align import substitution_matrices
R = sys.argv[1]; MINLEN = int(sys.argv[2]); accs = sys.argv[3:]
CODE = {}; bases = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i = 0
for a in bases:
    for b in bases:
        for c in bases: CODE[a+b+c] = aas[i]; i += 1
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
def orfs(g, mode, minlen):
    for c, s in g.items():
        L = len(s)
        for strand, seq in ((1, s), (-1, s.translate(COMP)[::-1])):
            for fr in range(3):
                prot = []; start = fr
                for p in range(fr, L - 2, 3):
                    cod = seq[p:p+3]; aa = CODE.get(cod, 'X')
                    if mode == 'tagsense' and cod == 'TAG': aa = 'X'
                    if aa == '*':
                        if len(prot) >= minlen: yield c, strand, start, p, ''.join(prot), True
                        prot = []; start = p + 3
                    else: prot.append(aa)
                if len(prot) >= minlen: yield c, strand, start, start + 3 * len(prot), ''.join(prot), False
def coords(strand, start, end, L):
    return (start + 1, end) if strand == 1 else (L - end + 1, L - start)
alpha = Alphabet.amino()
with HMMFile(f'{R}/data/rf12_sprot.hmm') as hf: hmms = list(hf)
ref = {k.split('|')[0]: v for k, v in read_fasta(f'{R}/data/rf_queries.faa').items()}
ec1, ec2 = ref['Ec_prfA'], ref['Ec_prfB']
i_pxt = ec1.find('PAT', 150) + 1; i_spf = ec2.find('SPF') + 1; i_ggq = ec1.find('GGQ') + 1
print('Ec RF1 PxT at', i_pxt, ec1[i_pxt-1:i_pxt+2], '| Ec RF2 SPF at', i_spf, ec2[i_spf-1:i_spf+2], '| GGQ at', i_ggq, file=sys.stderr)
aligner = Align.PairwiseAligner(); aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -11; aligner.extend_gap_score = -1; aligner.mode = 'local'
def map_motif(orf, refseq, start):
    aln = aligner.align(refseq, orf.replace('X', 'A'))[0]
    m = {}
    for (rs, re_), (qs, qe) in zip(*aln.aligned):
        for k in range(re_ - rs): m[rs + k + 1] = orf[qs + k]
    return ''.join(m.get(p, '-') for p in (start, start + 1, start + 2)), aln.score
out = open(f'{R}/out/rf_hmm.tsv', 'w')
out.write('genome\tmode\tcontig\tstrand\tstart\tend\taa_len\tends_in_stop\tRF1_score\tRF2_score\tRF1_E\tRF2_E\tcall\tEcRF1_PxT\tEcRF2_SPF\tGGQ\tsw_vs_EcRF1\tsw_vs_EcRF2\tX_count\n')
fa = open(f'{R}/out/rf_candidates.faa', 'w')
for acc in accs:
    g = read_fasta(f'{R}/data/{acc}.fna')
    for mode in ('std', 'tagsense'):
        recs = []; seqs = []
        for i, rec in enumerate(orfs(g, mode, MINLEN)):
            recs.append(rec)
            seqs.append(TextSequence(name=f'{acc}|{mode}|{i}'.encode(), sequence=rec[4]).digitize(alpha))
        sb = DigitalSequenceBlock(alpha, seqs)
        score = collections.defaultdict(dict)
        for th in pyhmmer.hmmsearch(hmms, sb, cpus=3, E=1e-3, Z=1000):
            hname = th.query.name; hname = hname.decode() if isinstance(hname, bytes) else hname
            for h in th:
                n = h.name.decode() if isinstance(h.name, bytes) else h.name
                score[n][hname] = (round(h.score, 1), h.evalue)
        for n, sc in sorted(score.items(), key=lambda x: int(x[0].split('|')[2])):
            i = int(n.split('|')[2]); c, strand, st, en, prot, stp = recs[i]
            s1, e1 = coords(strand, st, en, len(g[c]))
            a, ea = sc.get('RF1_sprot', (0, 1)); b, eb = sc.get('RF2_sprot', (0, 1))
            call = 'RF1' if a - b >= 50 else 'RF2' if b - a >= 50 else 'ambiguous'
            pxt, sw1 = map_motif(prot, ec1, i_pxt); spf, sw2 = map_motif(prot, ec2, i_spf); ggq, _ = map_motif(prot, ec1, i_ggq)
            out.write(f'{acc}\t{mode}\t{c}\t{strand}\t{s1}\t{e1}\t{len(prot)}\t{stp}\t{a}\t{b}\t{ea:.1e}\t{eb:.1e}\t{call}\t{pxt}\t{spf}\t{ggq}\t{sw1:.0f}\t{sw2:.0f}\t{prot.count("X")}\n')
            fa.write(f'>{acc}|{mode}|{c}|{s1}-{e1}|{strand}\n{prot}\n')
    print('done', acc, flush=True)
out.close(); fa.close()
