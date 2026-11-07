"""Calibrated version of ortholog_codons: for EVERY codon of the target genome (not just the test codon), the
distribution of the orthologous residue in close relatives. Then the test codon's profile is compared with each
sense codon's profile (which amino acid's codons does it look like?) and with the background (all codons).
Only orthologs with identity >= MINID; one vote per site from the most similar relative.
usage: python src/ortholog_profiles.py <target_acc> <codon> <table> <rel_acc> ..."""
import math
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import CODE
from ortholog_codons import genes

AAS = "ACDEFGHIKLMNPQRSTVWY"
MINID = float(__import__("os").environ.get("MINID", 60))


def main(target, test, table, rels):
    tg = genes(target, table)
    tprot = {k: "".join("X" if nt[i:i + 3] == test else CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 3, 3))
             for k, nt in tg.items()}
    best = {}
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.faa"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in tprot.items()))
        for rel in rels:
            rg = genes(rel, 11)
            db = Path(td) / f"{rel}.faa"
            db.write_text("".join(f">{k}\n{''.join(CODE.get(nt[i:i+3], 'X') for i in range(0, len(nt) - 3, 3))}\n"
                                  for k, nt in rg.items()))
            subprocess.run(["diamond", "makedb", "--in", str(db), "-d", str(db), "--quiet"], check=True)
            h = Path(td) / f"{rel}.tsv"
            subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(db) + ".dmnd", "-o", str(h), "--threads", "2",
                            "--quiet", "--evalue", "1e-10", "--max-target-seqs", "1", "--max-hsps", "1",
                            "--outfmt", "6", "qseqid", "pident", "qstart", "qseq", "sseq"], check=True)
            for line in open(h):
                qid, pid, qs, qseq, sseq = line.rstrip("\n").split("\t")
                pid = float(pid)
                if pid < MINID:
                    continue
                qi = int(qs) - 1
                for a, b in zip(qseq, sseq):
                    if a != "-":
                        if b in AAS:
                            k = (qid, qi)
                            if k not in best or pid > best[k][0]:
                                best[k] = (pid, b)
                        qi += 1
    prof = defaultdict(Counter)
    for (qid, qi), (pid, b) in best.items():
        prof[tg[qid][3 * qi:3 * qi + 3]][b] += 1
    bg = Counter()
    for c, v in prof.items():
        if c != test and CODE.get(c) not in ("*", None):
            bg.update(v)
    nb = sum(bg.values())
    t = prof[test]
    nt = sum(t.values())
    print(f"{target}: test codon {test}, orthologs >= {MINID}% identity; {nt} test sites, {nb} background sites")
    print("  residue  test%   background%  enrichment")
    for r, k in t.most_common(10):
        print(f"    {r}     {100*k/nt:5.1f}   {100*bg[r]/nb:5.1f}       {k/nt/(bg[r]/nb):.2f}")
    # which codon class does the test profile resemble? log-likelihood of test sites under each amino acid's profile
    aa_prof = defaultdict(Counter)
    for c, v in prof.items():
        a = CODE.get(c)
        if a not in ("*", None) and c != test:
            aa_prof[a].update(v)
    ll = []
    for a, v in aa_prof.items():
        tot = sum(v.values()) + 20
        ll.append((sum(k * math.log((v[r] + 1) / tot) for r, k in t.items()), a))
    tot = nb + 20
    ll.append((sum(k * math.log((bg[r] + 1) / tot) for r, k in t.items()), "background(any codon)"))
    ll.sort(reverse=True)
    print("  log-likelihood of the test sites' orthologous residues, relative to best:")
    print("   ", ", ".join(f"{a} {v - ll[0][0]:+.0f}" for v, a in ll[:8]))
    # mixture: P(r | test) = sum_A pi_A P(r | A) + pi_bg P(r | background); EM for pi
    comps = {a: {r: (aa_prof[a][r] + 1) / (sum(aa_prof[a].values()) + 20) for r in AAS} for a in aa_prof}
    comps["bg"] = {r: (bg[r] + 1) / (nb + 20) for r in AAS}
    pi = {a: 1 / len(comps) for a in comps}
    for _ in range(500):
        acc = {a: 0.0 for a in comps}
        for r, k in t.items():
            z = sum(pi[a] * comps[a][r] for a in comps)
            for a in comps:
                acc[a] += k * pi[a] * comps[a][r] / z
        pi = {a: acc[a] / nt for a in comps}
    llm = sum(k * math.log(sum(pi[a] * comps[a][r] for a in comps)) for r, k in t.items())
    print("  mixture weights (EM):", ", ".join(f"{a} {v:.2f}" for a, v in sorted(pi.items(), key=lambda x: -x[1]) if v >= 0.02),
          f"  [logL {llm - ll[0][0]:+.0f} vs best single]")
    for c in ("GCC", "GCG", "CAG", "AAG", "TGG"):
        v = prof[c]
        n = sum(v.values())
        if n:
            print(f"  calibration {c} ({CODE[c]}): n={n}, orthologous residue {', '.join(f'{r} {100*k/n:.0f}%' for r, k in v.most_common(3))}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4:])
