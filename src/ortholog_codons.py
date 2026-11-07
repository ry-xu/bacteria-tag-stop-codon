"""For each in-frame test codon (e.g. TAG) in a target genome, which codon do close relatives have at the
orthologous position?

Target genes are called with a table that reads the test codon as sense (15 for TAG), translated with the test
codon as 'X', and aligned (DIAMOND) to the proteins of each relative genome (table 11). For every target gene the
best hit per relative is kept (identity >= MINID), and at each aligned 'X' the relative's codon is recorded.
A codon with one meaning A should sit where relatives have A codons; a codon that just arose by point mutation
sits where relatives have its one-step neighbours (CAG, AAG, GAG, TGG, TTG, TCG, TAT/TAC, TAA/TGA).

usage: python src/ortholog_codons.py <target_acc> <codon> <table> <rel_acc> [<rel_acc> ...]
       -> results/ortho/<target>_<codon>.tsv (one row per site x relative) and a summary on stdout
"""
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import CODE, ROOT, fetch, read_fasta

MINID = 50.0


def genes(acc, table):
    contigs = dict(read_fasta(fetch(acc)))
    total = sum(len(s) for s in contigs.values())
    finder = pyrodigal.GeneFinder(meta=total < 100000)
    if total >= 100000:
        finder.train(*[s.encode() for s in contigs.values()], translation_table=table)
    out = {}
    for c, s in contigs.items():
        for g in finder.find_genes(s.encode()):
            out[f"{c}|{g.begin}|{g.end}|{g.strand}"] = g.sequence()
    return out


def main(target, test, table, rels):
    tg = genes(target, table)
    tprot = {k: "".join("X" if nt[i:i + 3] == test else CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 3, 3))
             for k, nt in tg.items()}
    tprot = {k: v for k, v in tprot.items() if "X" in v}
    rows = []
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.faa"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in tprot.items()))
        for rel in rels:
            rel, rt = (rel.split(":")[0], int(rel.split(":")[1])) if ":" in rel else (rel, 11)
            rg = genes(rel, rt)
            db = Path(td) / f"{rel}.faa"
            db.write_text("".join(f">{k}\n{''.join(CODE.get(nt[i:i+3], 'X') for i in range(0, len(nt) - 3, 3))}\n"
                                  for k, nt in rg.items()))
            subprocess.run(["diamond", "makedb", "--in", str(db), "-d", str(db), "--quiet"], check=True)
            h = Path(td) / f"{rel}.tsv"
            subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(db) + ".dmnd", "-o", str(h), "--threads", "2",
                            "--quiet", "--evalue", "1e-10", "--max-target-seqs", "1", "--max-hsps", "1",
                            "--outfmt", "6", "qseqid", "sseqid", "pident", "qstart", "sstart", "qseq", "sseq"], check=True)
            for line in open(h):
                qid, sid, pid, qs, ss, qseq, sseq = line.rstrip("\n").split("\t")
                if float(pid) < MINID:
                    continue
                qi, si = int(qs) - 1, int(ss) - 1
                for a, b in zip(qseq, sseq):
                    if a == "X" and b != "-":
                        rcod = rg[sid][3 * si:3 * si + 3]
                        rows.append((qid, qi, rel, float(pid), rcod, CODE.get(rcod, "?"), b))
                    if a != "-":
                        qi += 1
                    if b != "-":
                        si += 1
    out = ROOT / "results" / "ortho"
    out.mkdir(parents=True, exist_ok=True)
    with open(out / f"{target}_{test}.tsv", "w") as f:
        f.write("gene\taa_pos\trelative\tpident\trel_codon\trel_aa\trel_residue\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    n_sites = sum(v.count("X") for v in tprot.values())
    print(f"{target}: {n_sites} in-frame {test} in {len(tprot)} genes; {len(rows)} site x relative pairs")
    # one vote per site: the codon in the most similar relative
    best = {}
    for qid, qi, rel, pid, rcod, raa, b in rows:
        k = (qid, qi)
        if k not in best or pid > best[k][0]:
            best[k] = (pid, rcod, raa)
    aa = Counter(v[2] for v in best.values())
    cod = Counter(v[1] for v in best.values())
    print(f"sites with an orthologous codon: {len(best)}")
    print("  relative's amino acid:", aa.most_common())
    print("  relative's codon     :", cod.most_common(20))
    # sites where all relatives with a hit agree on the amino acid
    per = defaultdict(list)
    for qid, qi, rel, pid, rcod, raa, b in rows:
        per[(qid, qi)].append(raa)
    agree = Counter(v[0] for v in per.values() if len(v) >= 3 and len(set(v)) == 1)
    print("  sites where >=3 relatives all have the same amino acid:", agree.most_common())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4:])
