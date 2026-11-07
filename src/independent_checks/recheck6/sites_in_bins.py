"""Do the two independent assemblies (bin.237, bin.296) carry the candidate's in-frame TAGs? blastn the +-75 nt
window of every candidate TAG site (out/crossmap_windows.fna, written by crossmap_check.py) against each bin and read
the codon at the site from the best hit (>= 100 nt, >= 95% id).  usage: sites_in_bins.py"""
import subprocess
from collections import Counter
from pathlib import Path
D = Path(__file__).resolve().parent.parent
q = D / "out" / "crossmap_windows.fna"
for b in ("QTP_TA20190011_bin.237", "QTP_TA20190013_bin.296"):
    out = subprocess.run(["blastn", "-query", str(q), "-db", str(D / "out" / "blastdb" / b), "-evalue", "1e-20",
                          "-max_target_seqs", "1", "-max_hsps", "1", "-outfmt", "6 qseqid pident length qstart qend qseq sseq"],
                         capture_output=True, text=True, check=True).stdout
    cod = {}
    for line in out.splitlines():
        qid, pid, ln, qs, qe, qa, sa = line.split("\t")
        off = int(qid.split("|")[3])
        if int(ln) < 100 or float(pid) < 95 or qid in cod or not (int(qs) - 1 <= off and off + 3 <= int(qe)):
            continue
        qi, c = int(qs) - 2, ""
        for a, s in zip(qa, sa):
            if a != "-":
                qi += 1
            if off <= qi < off + 3 and a != "-":
                c += s
        cod[qid] = c
    n = sum(1 for l in open(q) if l[0] == ">")
    print(f"{b}: candidate TAG sites {n}; found in bin {len(cod)}; codon there: {Counter(cod.values()).most_common(6)}")
