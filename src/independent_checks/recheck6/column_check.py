"""Are the TAG sites real conserved-Trp columns? For each TAG site of a genome (from out/stops/<label>.sites.tsv),
translate +-60 codons around it (code 11, the TAG masked as X, other TAGs X, cut at TAA/TGA), DIAMOND blastp vs
Swiss-Prot (top 50, e<=1e-5), and read the residue each hit has at the masked column. Calibration: the same for
200 random TGG codons (masked as X) inside the genome's own table-11 genes.
Also merges the residue that each relative proteome puts at the site (out/stops/<genome>_by*.sites.tsv).
usage: column_check.py <genome> <sites_label> -> out/columns/<genome>.tsv + summary"""
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

D = Path(__file__).resolve().parent.parent
SPROT = D.parent / "data" / "sprot.dmnd"
COMP = str.maketrans("ACGTN", "TGCAN")
CODE = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODE[a + b + c] = aas[i]; i += 1


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


def window(dna, c, p0, st, flank=60):
    s = dna[c]
    cods = []
    for j in range(-flank, flank + 1):
        q = p0 + 3 * j if st == 1 else p0 - 3 * j
        if q < 0 or q + 3 > len(s):
            cods.append(None); continue
        cod = s[q:q + 3] if st == 1 else rc(s[q:q + 3])
        cods.append(cod)
    # left part: walk back from the centre until a non-TAG stop / edge
    aa = []
    centre = flank
    left = centre
    while left - 1 >= 0 and cods[left - 1] and (CODE.get(cods[left - 1], "X") != "*" or cods[left - 1] == "TAG"):
        left -= 1
    right = centre
    while right + 1 < len(cods) and cods[right + 1] and (CODE.get(cods[right + 1], "X") != "*" or cods[right + 1] == "TAG"):
        right += 1
    prot = "".join("X" if (k == centre or cods[k] == "TAG") else CODE.get(cods[k], "X") for k in range(left, right + 1))
    return prot, centre - left


def main(genome, label):
    dna = dict(read_fasta(D / "data" / f"{genome}.fna"))
    sites = []
    for line in list(open(D / "out" / "stops" / f"{label}.sites.tsv"))[1:]:
        x = line.rstrip("\n").split("\t")
        if x[3] == "TAG":
            sites.append(("TAG", x[0], int(x[1]), int(x[2])))
    # relatives' residues per site
    rel = defaultdict(dict)
    for q, lab in (("CAKU1", f"{genome}_byQ1"), ("CAKU2", f"{genome}_byQ2"), ("relB", f"{genome}_byRelB")):
        p = D / "out" / "stops" / f"{lab}.sites.tsv"
        if p.exists():
            for line in list(open(p))[1:]:
                x = line.rstrip("\n").split("\t")
                rel[(x[0], int(x[1]))][q] = x[4]
    # TGG calibration: random internal TGG codons of table-11 genes
    random.seed(1)
    genes = [l.rstrip("\n").split("\t") for l in list(open(D / "out" / "genes" / f"{genome}.tsv"))[1:]]
    cal = []
    while len(cal) < 200:
        g = random.choice(genes)
        c, b, e, st = g[1], int(g[2]), int(g[3]), int(g[4])
        if g[6] != "00" or e - b < 600:
            continue
        nt = dna[c][b - 1:e] if st == 1 else rc(dna[c][b - 1:e])
        tgg = [k for k in range(60 * 3, len(nt) - 60 * 3, 3) if nt[k:k + 3] == "TGG"]
        if not tgg:
            continue
        k = random.choice(tgg)
        p0 = b - 1 + k if st == 1 else e - 3 - k
        cal.append(("TGG", c, p0, st))
    out = D / "out" / "columns"
    out.mkdir(exist_ok=True)
    faa = out / f"{genome}.query.faa"
    pos = {}
    with open(faa, "w") as fh:
        for i, (kind, c, p0, st) in enumerate(sites + cal):
            prot, m = window(dna, c, p0, st)
            if len(prot) < 40:
                continue
            sid = f"{kind}|{c}|{p0}|{st}"
            pos[sid] = m
            fh.write(f">{sid}\n{prot}\n")
    hits = out / f"{genome}.sprot.tsv"
    if not hits.exists():
        subprocess.run(["diamond", "blastp", "-q", str(faa), "-d", str(SPROT), "-o", str(hits), "--more-sensitive",
                        "--max-target-seqs", "50", "--evalue", "1e-5", "--threads", "4", "--quiet", "--outfmt", "6",
                        "qseqid", "sseqid", "pident", "qstart", "qend", "evalue", "qseq_gapped", "sseq_gapped"], check=True)
    res = defaultdict(Counter)
    for line in open(hits):
        q, s, pid, qs, qe, ev, qa, sa = line.rstrip("\n").split("\t")
        qs, qe = int(qs), int(qe)
        m = pos[q]
        if not (qs - 1 <= m <= qe - 1):
            continue
        qi = qs - 2
        for a, b in zip(qa, sa):
            if a != "-":
                qi += 1
            if qi == m and a != "-":
                res[q][b] += 1
                break
    rows = []
    for sid in pos:
        r = res.get(sid, Counter())
        n = sum(r.values())
        rows.append((sid, n, r))
    with open(out / f"{genome}.tsv", "w") as fh:
        fh.write("site\tn_sprot\tW_frac\tsprot_residues\tCAKU1\tCAKU2\trelB\n")
        for sid, n, r in rows:
            k = sid.split("|")
            rr = rel.get((k[1], int(k[2])), {})
            fh.write(f"{sid}\t{n}\t{(r['W'] / n if n else float('nan')):.2f}\t{' '.join(f'{a}{c}' for a, c in r.most_common(5))}"
                     f"\t{rr.get('CAKU1', '.')}\t{rr.get('CAKU2', '.')}\t{rr.get('relB', '.')}\n")
    for kind in ("TAG", "TGG"):
        rs = [(n, r) for sid, n, r in rows if sid.startswith(kind) and n >= 5]
        maj = Counter(r.most_common(1)[0][0] for n, r in rs)
        wcons = sum(r["W"] / n >= 0.8 for n, r in rs)
        wsome = sum(0.2 <= r["W"] / n < 0.8 for n, r in rs)
        meanw = sum(r["W"] / n for n, r in rs) / max(1, len(rs))
        print(f"{genome} {kind}: {len(rs)} sites with >= 5 Swiss-Prot hits at the column; majority residue {maj.most_common(8)}; "
              f"W in >= 80% of hits: {wcons} ({wcons / max(1, len(rs)):.0%}); W in 20-80%: {wsome}; mean W fraction {meanw:.2f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
