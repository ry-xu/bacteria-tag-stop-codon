"""Count stop codons that sit INSIDE high-identity tblastn alignments (relative's protein vs a genome's raw DNA).

For each target genome, combine the tblastn tables of the chosen query proteomes (out/tbn/<Q>__<T>.tsv). An HSP
qualifies if pident >= MINPID and length >= 60 aa. A '*' column qualifies if it has >= FLANK aligned columns on
each side inside the HSP and the 30-column windows on each side are each >= WINID identical. Each stop site
(contig, codon start, strand) is counted once (best-scoring HSP kept). The codon is re-read from the DNA (strand
aware) and checked to be a stop. Also records the query residue at the stop column (for meaning) and, as a
calibration, the query residue opposite every TGG / CAG / TAC / TAT codon inside the same qualifying windows.

usage: stops_in_hsps.py <label> <target> <Q1,Q2,...> [minpid=50]  -> out/stops/<label>.sites.tsv, prints a summary
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

D = Path(__file__).resolve().parent.parent
FLANK, WIN, WINID = 20, 30, 0.40
CAL = {"TGG": "W", "CAG": "Q", "TAC": "Y", "TAT": "Y", "GAG": "E", "AAG": "K"}
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(p):
    name, seq = None, []
    for line in open(p):
        if line[0] == ">":
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip().upper())
    if name:
        yield name, "".join(seq)


def main(label, target, queries, minpid=50.0):
    dna = dict(read_fasta(D / "data" / f"{target}.fna"))
    sites = {}
    cal = defaultdict(Counter)       # codon -> query residues opposite it (deduplicated by site)
    cal_seen = {}
    aligned_codons = set()
    bad = 0
    for q in queries:
        for line in open(D / "out" / "tbn" / f"{q}__{target}.tsv"):
            f = line.rstrip("\n").split("\t")
            qid, sid, pid, ln = f[0], f[1], float(f[2]), int(f[3])
            qs, qe, qlen, ss, se, sfr, bits = int(f[4]), int(f[5]), int(f[6]), int(f[7]), int(f[8]), int(f[9]), float(f[11])
            qa, sa = f[12], f[13]
            if pid < minpid or ln < 60:
                continue
            strand = 1 if sfr > 0 else -1
            n = len(qa)
            ident = [1 if (a == b and a != "-") else 0 for a, b in zip(qa, sa)]
            cum = [0]
            for x in ident:
                cum.append(cum[-1] + x)
            k = -1          # subject aa index
            qp = qs - 1     # query position (0-based) of last non-gap query residue
            for col in range(n):
                if qa[col] != "-":
                    qp += 1
                if sa[col] == "-":
                    continue
                k += 1
                if col < FLANK or n - 1 - col < FLANK:
                    continue
                lo, hi = max(0, col - WIN), min(n, col + 1 + WIN)
                lid = (cum[col] - cum[lo]) / (col - lo)
                rid = (cum[hi] - cum[col + 1]) / (hi - col - 1)
                if lid < WINID or rid < WINID:
                    continue
                if strand == 1:
                    p0 = ss - 1 + 3 * k
                    cod = dna[sid][p0:p0 + 3]
                else:
                    p0 = ss - 3 - 3 * k          # 0-based start (forward coords) of the codon on the minus strand
                    cod = rc(dna[sid][p0:p0 + 3])
                key = (sid, p0, strand)
                aligned_codons.add(key)
                if sa[col] == "*":
                    if cod not in ("TAA", "TAG", "TGA"):
                        bad += 1
                        continue
                    rec = (bits, cod, qa[col], qid, qp + 1, qlen, pid, round(lid, 2), round(rid, 2))
                    if key not in sites or bits > sites[key][0]:
                        sites[key] = rec
                elif cod in CAL and key not in cal_seen:
                    cal_seen[key] = 1
                    cal[cod][qa[col]] += 1
    out = D / "out" / "stops"
    out.mkdir(exist_ok=True)
    with open(out / f"{label}.sites.tsv", "w") as fh:
        fh.write("contig\tpos0\tstrand\tcodon\tquery_res\tquery\tqpos\tqlen\thsp_pid\tleft_id\tright_id\tbits\n")
        for (c, p0, st), r in sorted(sites.items()):
            fh.write(f"{c}\t{p0}\t{st}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]}\t{r[8]}\t{r[0]}\n")
    cnt = Counter(r[1] for r in sites.values())
    res = Counter(r[2] for r in sites.values() if r[1] == "TAG")
    calline = "; ".join(f"{c}->{a} {cal[c][a] / max(1, sum(cal[c].values())):.0%} (n={sum(cal[c].values())})"
                        for c, a in CAL.items())
    print(f"{label}\tqueries={','.join(queries)}\tminpid={minpid}\taligned_codons={len(aligned_codons)}\t"
          f"TAG={cnt['TAG']}\tTAA={cnt['TAA']}\tTGA={cnt['TGA']}\tnonstop_bug={bad}")
    print(f"  query residue at TAG: {res.most_common()}")
    print(f"  calibration: {calline}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3].split(","), float(sys.argv[4]) if len(sys.argv) > 4 else 50.0)
