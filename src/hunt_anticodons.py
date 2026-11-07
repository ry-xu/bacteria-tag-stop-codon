"""Tally anticodons in reads that match the clade's tRNAs (output of results/trna_hunt/run_hunt.sh).

For each SAM hit (read aligned to a clade tRNA, minimap2 -x sr), require >= MINLEN aligned bases at >= MINID identity
outside the anticodon, then read the 3 bases aligned to the reference anticodon. A tRNA whose body matches one of the
organism's tRNAs but whose anticodon is CTA would be a candidate UAG reader derived from that tRNA.

usage: python src/hunt_anticodons.py results/trna_hunt/part_*.sam
"""
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MINLEN, MINID = 55, float(__import__("os").environ.get("MINID", 0.95))
COMP = str.maketrans("ACGTN", "TGCAN")


def main(paths):
    meta = {r["name"]: r for r in csv.DictReader(open(ROOT / "results" / "trna_hunt" / "clade_trnas.tsv"), delimiter="\t")}
    ref = {}
    name = None
    for l in open(ROOT / "results" / "trna_hunt" / "clade_trnas.fa"):
        if l.startswith(">"):
            name = l[1:].strip()
            ref[name] = ""
        else:
            ref[name] += l.strip()
    tally = defaultdict(Counter)
    examples = {}
    n = kept = 0
    for p in paths:
        for line in open(p):
            if line.startswith("@"):
                continue
            f = line.split("\t")
            if len(f) < 11:
                continue  # partial line (file still being written)
            flag, rname, pos, cigar, seq = int(f[1]), f[2], int(f[3]), f[5], f[9]
            if flag & 4 or rname not in meta or seq == "*":
                continue
            n += 1
            a0, a1 = int(meta[rname]["ac_start"]) - 1, int(meta[rname]["ac_end"])  # 0-based anticodon span on ref
            r, q = pos - 1, 0
            match = mism = 0
            ac = {}
            for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                ln = int(ln)
                if op in "M=X":
                    if q + ln > len(seq) or r + ln > len(ref[rname]):
                        break
                    for k in range(ln):
                        rb, qb = ref[rname][r + k], seq[q + k]
                        if a0 <= r + k < a1:
                            ac[r + k] = qb
                        elif rb == qb:
                            match += 1
                        else:
                            mism += 1
                    r += ln
                    q += ln
                elif op in "IS":
                    q += ln
                elif op in "DN":
                    for k in range(ln):
                        if a0 <= r + k < a1:
                            ac[r + k] = "-"
                    mism += ln
                    r += ln
            if match + mism < MINLEN or match / (match + mism) < MINID or len(ac) < 3:
                continue
            kept += 1
            triplet = "".join(ac[i] for i in range(a0, a1))
            tally[rname][triplet] += 1
            if triplet != meta[rname]["anticodon"] and (rname, triplet) not in examples:
                examples[rname, triplet] = seq
    print(f"hits {n}, kept {kept} (>= {MINLEN} nt at >= {MINID:.0%} outside the anticodon)")
    for rname in sorted(tally, key=lambda k: -sum(tally[k].values())):
        c = tally[rname]
        own = meta[rname]["anticodon"]
        other = {k: v for k, v in c.items() if k != own}
        flag = "  <== CTA" if "CTA" in other else ""
        print(f"  {rname:45s} own {own}: {c[own]:5d}   others: {dict(Counter(other).most_common(4))}{flag}")
    for (rname, t), s in examples.items():
        if t == "CTA":
            print("CTA example read for", rname, ":", s)


if __name__ == "__main__":
    main(sys.argv[1:])
