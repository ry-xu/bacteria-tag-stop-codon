"""Meaning of a test codon by conservation stratum (cf. codon-hunt's lesson: variable columns report where a codon
came from, conserved columns report what it means). Homologs = top 25 Swiss-Prot hits (DIAMOND). For each column
with >= MINH homologs, agreement = share of the majority residue. For each stratum, P(majority residue | meaning A)
is estimated from the genome's own codons for A, and the test codon is scored under every A.
usage: python src/codon_strata.py <acc> <table> <codon>"""
import math
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import CODE, DB, fetch, read_fasta

AAS = "ACDEFGHIKLMNPQRSTVWY"
MINH = 8
STRATA = [(0.0, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]
acc, table, test = sys.argv[1], int(sys.argv[2]), sys.argv[3]
contigs = dict(read_fasta(fetch(acc)))
finder = pyrodigal.GeneFinder()
finder.train(*[s.encode() for s in contigs.values()], translation_table=table)
cds = {}
for c, s in contigs.items():
    for g in finder.find_genes(s.encode()):
        cds[f"{c}|{g.begin}|{g.end}|{g.strand}"] = g.sequence()
with tempfile.TemporaryDirectory() as td:
    q = Path(td) / "q.faa"
    q.write_text("".join(f">{k}\n" + "".join("X" if nt[i:i+3] == test else CODE.get(nt[i:i+3], "X")
                                             for i in range(0, len(nt) - 3, 3)) + "\n" for k, nt in cds.items()))
    h = Path(td) / "h.tsv"
    subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(DB), "-o", str(h), "--threads", "4", "--quiet",
                    "--evalue", "1e-10", "--max-target-seqs", "25", "--max-hsps", "1",
                    "--outfmt", "6", "qseqid", "qstart", "qseq", "sseq"], check=True)
    hits = defaultdict(list)
    for line in open(h):
        qid, qs, qseq, sseq = line.rstrip("\n").split("\t")
        hits[qid].append((int(qs), qseq, sseq))
obs = defaultdict(list)  # stratum -> [(codon, majority residue)]
for qid, hl in hits.items():
    nt = cds[qid]
    cols = defaultdict(list)
    for qs, qseq, sseq in hl:
        qi = qs - 1
        for a, b in zip(qseq, sseq):
            if a != "-":
                if b in AAS:
                    cols[qi].append(b)
                qi += 1
    for qi, res in cols.items():
        if len(res) < MINH:
            continue
        top, k = Counter(res).most_common(1)[0]
        agr = k / len(res)
        for si, (lo, hi) in enumerate(STRATA):
            if lo <= agr < hi:
                obs[si].append((nt[3 * qi:3 * qi + 3], top))
print(f"{acc}  test codon {test}  (homologs >= {MINH}; majority residue per column)")
for si, (lo, hi) in enumerate(STRATA):
    pairs = obs[si]
    m = {}
    for A in AAS:
        c = Counter(r for cod, r in pairs if CODE.get(cod) == A and cod != test)
        tot = sum(c.values()) + len(AAS)
        m[A] = {r: (c[r] + 1) / tot for r in AAS}
    t = [r for cod, r in pairs if cod == test]
    all_n = sum(1 for cod, r in pairs if CODE.get(cod) not in ("*", None))
    ll = sorted(((sum(math.log(m[A][r]) for r in t), A) for A in AAS), reverse=True)
    rank = ", ".join(f"{A}{v - ll[0][0]:+.1f}" for v, A in ll[:5])
    print(f"  agreement {lo:.1f}-{min(hi,1):.1f}: test n={len(t):4d} ({len(t)/max(all_n,1)*1000:.2f} per 1000 codons)  "
          f"best: {rank}\n      majority residues: {Counter(t).most_common(8)}")
