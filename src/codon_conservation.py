"""Per-codon share of occurrences at strictly conserved alignment columns (all Swiss-Prot top hits agree), and which
residue the homologs hold there. A codon with one clean meaning sits at conserved sites about as often as the other
codons for that amino acid; an ambiguously decoded codon is pushed out of them.
usage: python src/codon_conservation.py <acc> <table> <codon>"""
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import CODE, DB, fetch, read_fasta

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
                    "--evalue", "1e-10", "--max-target-seqs", "10", "--max-hsps", "1",
                    "--outfmt", "6", "qseqid", "qstart", "qseq", "sseq"], check=True)
    hits = defaultdict(list)
    for line in open(h):
        qid, qs, qseq, sseq = line.rstrip("\n").split("\t")
        hits[qid].append((int(qs), qseq, sseq))
tot, cons, cons_res = Counter(), Counter(), defaultdict(Counter)
for qid, hl in hits.items():
    nt = cds[qid]
    cols = defaultdict(list)
    for qs, qseq, sseq in hl:
        qi = qs - 1
        for a, b in zip(qseq, sseq):
            if a != "-":
                if b != "-":
                    cols[qi].append(b)
                qi += 1
    for qi, res in cols.items():
        if len(res) < 5:
            continue
        cod = nt[3 * qi:3 * qi + 3]
        tot[cod] += 1
        if len(set(res)) == 1:
            cons[cod] += 1
            cons_res[cod][res[0]] += 1
rows = sorted(((cons[c] / tot[c], c, CODE.get(c, "?"), tot[c], cons[c]) for c in tot if tot[c] >= 20))
print(f"{acc}: share of codon occurrences at strictly conserved columns (>=5 homologs, all agree)")
for frac, c, aa, n, k in rows:
    mark = "  <==" if c == test else ""
    top = ",".join(f"{r}{v}" for r, v in cons_res[c].most_common(4))
    print(f"  {c} {aa if c != test else 'X'}  {frac:.3f}  n={n:5d} conserved={k:4d}  [{top}]{mark}")
