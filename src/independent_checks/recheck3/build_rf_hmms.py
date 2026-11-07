"""Build full-length RF1 (prfA) and RF2 (prfB) profile HMMs locally (no download) from Swiss-Prot entries named exactly
'Peptide chain release factor 1|RF1' / 'Peptide chain release factor 2|RF2'. Seed = single-sequence HMM of E. coli RF1/RF2
(from the Swiss-Prot set itself), then 3 rounds of hmmalign-all -> hmmbuild. Also writes the RF reference fasta and a
self-check: every Swiss-Prot RF1/RF2 is scored with both HMMs (discrimination margin).
Usage: uv run --no-project --with pyhmmer python -I build_rf_hmms.py <recheck3_dir> <uniprot_sprot.fasta.gz>"""
import sys, gzip, re
import pyhmmer
from pyhmmer.easel import TextSequence, Alphabet, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background
R, SP = sys.argv[1], sys.argv[2]
alpha = Alphabet.amino(); bg = Background(alpha)
sets = {'RF1': [], 'RF2': []}
name = None; desc = None; buf = []
def flush():
    if name is None: return
    d = desc.split(' OS=')[0]
    if d in ('Peptide chain release factor 1', 'Peptide chain release factor RF1'): sets['RF1'].append((name, ''.join(buf)))
    elif d in ('Peptide chain release factor 2', 'Peptide chain release factor RF2'): sets['RF2'].append((name, ''.join(buf)))
with gzip.open(SP, 'rt') as f:
    for line in f:
        if line.startswith('>'):
            flush(); h = line[1:].rstrip(); name = h.split()[0].split('|')[2]; desc = h.split(' ', 1)[1]; buf = []
        else: buf.append(line.strip())
    flush()
with open(f'{R}/data/sprot_rf12.faa', 'w') as o:
    for k, v in sets.items():
        for n, s in v: o.write(f'>{k}|{n}\n{s}\n')
hmms = {}
for k, v in sets.items():
    v = [(n, s) for n, s in v if 'X' not in s and 'U' not in s and 'B' not in s and 'Z' not in s]
    seqs = [TextSequence(name=n.encode(), sequence=s).digitize(alpha) for n, s in v]
    seed = [s for s, (n, _) in zip(seqs, v) if n == ('RF1_ECOLI' if k == 'RF1' else 'RF2_ECOLI')][0]
    builder = Builder(alpha)
    hmm, _, _ = builder.build(seed, bg)
    for it in range(3):
        msa = pyhmmer.hmmer.hmmalign(hmm, seqs, trim=True, digitize=True)
        msa.name = f'{k}_sprot'.encode()
        hmm, _, _ = Builder(alpha, symfrac=0.5).build_msa(msa, bg)
        hmm.name = f'{k}_sprot'.encode()
    hmms[k] = hmm
    print(k, 'n_seq', len(seqs), 'model length', hmm.M, file=sys.stderr)
with open(f'{R}/data/rf12_sprot.hmm', 'wb') as o:
    for k in ('RF1', 'RF2'): hmms[k].write(o)
# self-check
allseq = [TextSequence(name=f'{k}|{n}'.encode(), sequence=s).digitize(alpha) for k, v in sets.items() for n, s in v]
blk = DigitalSequenceBlock(alpha, allseq)
sc = {}
for th in pyhmmer.hmmsearch(list(hmms.values()), blk, cpus=3, E=10):
    q = th.query.name; q = q.decode() if isinstance(q, bytes) else q
    for h in th:
        n = h.name.decode() if isinstance(h.name, bytes) else h.name
        sc.setdefault(n, {})[q] = h.score
for k in ('RF1', 'RF2'):
    marg = sorted((s.get('RF1_sprot', 0) - s.get('RF2_sprot', 0)) for n, s in sc.items() if n.startswith(k))
    own = sorted(s.get(f'{k}_sprot', 0) for n, s in sc.items() if n.startswith(k))
    print(f'{k}: n={len(marg)} own-score min/median {own[0]:.0f}/{own[len(own)//2]:.0f}; RF1-RF2 margin min {marg[0]:.0f} 5% {marg[len(marg)//20]:.0f} median {marg[len(marg)//2]:.0f} max {marg[-1]:.0f}')
