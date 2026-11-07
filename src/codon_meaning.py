"""What does a read-through stop codon mean? Homolog-residue likelihood test on one genome.

Genes are called with a translation table that reads the stop as sense (15 for TAG, 4 for TGA), translated with the
test codon as 'X', and aligned to Swiss-Prot (DIAMOND, best hit per gene). At every aligned position we record
(codon, homolog residue). For each amino acid A, P(homolog residue | codon means A) is estimated from the genome's
own codons for A; the test codon's homolog residues are then scored under each A. Positions are also split by how
conserved the column is among the top hits (all hits agree vs not), since codon-hunt found that only conserved
sites report meaning reliably.

usage: python src/codon_meaning.py <acc> <TAG|TGA|TAA> [table]  -> prints ranking, writes results/meaning/<acc>_<codon>.json
"""
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import CODE, DB, ROOT, fetch, read_fasta

AAS = "ACDEFGHIKLMNPQRSTVWY"


def main(acc, codon, table):
    contigs = dict(read_fasta(fetch(acc)))
    finder = pyrodigal.GeneFinder()
    finder.train(*[s.encode() for s in contigs.values()], translation_table=table)
    cds = {}
    for c, s in contigs.items():
        for g in finder.find_genes(s.encode()):
            nt = g.sequence()
            cds[f"{c}|{g.begin}|{g.end}|{g.strand}"] = nt
    out = ROOT / "results" / "meaning"
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.faa"
        with open(q, "w") as f:
            for k, nt in cds.items():
                aa = "".join("X" if nt[i:i + 3] == codon else CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 3, 3))
                f.write(f">{k}\n{aa}\n")
        h = Path(td) / "h.tsv"
        subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(DB), "-o", str(h), "--threads", "4", "--quiet",
                        "--evalue", "1e-10", "--max-target-seqs", "10", "--max-hsps", "1",
                        "--outfmt", "6", "qseqid", "sseqid", "qstart", "sstart", "qseq", "sseq", "bitscore"], check=True)
        hits = defaultdict(list)
        for line in open(h):
            qid, sid, qs, ss, qseq, sseq, bs = line.rstrip("\n").split("\t")
            hits[qid].append((int(qs), qseq, sseq))
    # per query position: list of homolog residues over all hits; best hit first
    pairs = []  # (codon, best-hit residue, conserved?)
    for qid, hl in hits.items():
        nt = cds[qid]
        cols = defaultdict(list)
        for qs, qseq, sseq in hl:
            qi = qs - 1
            for a, b in zip(qseq, sseq):
                if a != "-":
                    if b != "-" and b in AAS:
                        cols[qi].append(b)
                    qi += 1
        for qi, res in cols.items():
            cod = nt[3 * qi:3 * qi + 3]
            if len(res) < 3:
                continue
            top = Counter(res).most_common(1)[0]
            pairs.append((cod, res[0], top[1] == len(res), top[0]))
    # model: P(homolog residue | meaning A), from standard codons for A
    def model(conserved_only):
        m = {}
        for A in AAS:
            c = Counter(r for cod, r, cons, top in pairs if CODE.get(cod) == A and cod != codon
                        and (cons or not conserved_only))
            tot = sum(c.values()) + len(AAS)
            m[A] = {r: (c[r] + 1) / tot for r in AAS}
        return m
    res = {}
    for label, cons_only in (("all", False), ("conserved", True)):
        m = model(cons_only)
        obs = [r for cod, r, cons, top in pairs if cod == codon and (cons or not cons_only)]
        ll = {A: sum(math.log(m[A][r]) for r in obs) for A in AAS}
        best = sorted(ll.items(), key=lambda x: -x[1])
        res[label] = dict(n=len(obs), residues=Counter(obs).most_common(), ranking=[(A, round(v - best[0][1], 1)) for A, v in best[:6]])
        print(f"{acc} {codon} [{label}] n={len(obs)}  ranking (logL minus best): {res[label]['ranking']}")
        print("   homolog residues:", Counter(obs).most_common(12))
    (out / f"{acc}_{codon}.json").write_text(json.dumps(res))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 11)
