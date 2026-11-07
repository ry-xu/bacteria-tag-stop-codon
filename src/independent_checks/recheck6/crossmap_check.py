"""Could the non-TAG alleles at the candidate's TAG sites come from reads of OTHER species in the same animals
(e.g. its relative GCA_056207115.1 / bin.73, present in TA20190013 and TA20190008)? read_support.py keeps reads
with <= 3% mismatches over the aligned (unclipped) part, so a read from another species can pass where the two
genomes are locally >= ~97% identical. For every TAG site with reads: take the candidate's +-75 nt window, blastn
it against all 584 non-conspecific CNGB bins of the four animals, and record the best local identity and the codon
the other genome has at the site.
usage: crossmap_check.py -> out/crossmap.tsv + summary"""
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

D = Path(__file__).resolve().parent.parent
RS = D.parent / "results" / "read_support"
RUNS = ["CRR939056", "CRR939061", "CRR939060", "CRR939058"]
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


dna = dict(read_fasta(D / "data" / "GCA_056209055.1.fna"))
obs = defaultdict(dict)
strand = {}
for run in RUNS:
    for line in list(open(RS / f"GCA_056209055.1_{run}.sites.tsv"))[1:]:
        c, p, st, cods = line.rstrip("\n").split("\t")
        strand[(c, int(p))] = int(st)
        obs[(c, int(p))][run] = Counter({k: int(v) for k, v in (x.rsplit(":", 1) for x in cods.split())})
W = 75
q = D / "out" / "crossmap_windows.fna"
with open(q, "w") as fh:
    for (c, p), st in strand.items():
        lo, hi = max(0, p - W), min(len(dna[c]), p + 3 + W)
        w = dna[c][lo:hi]
        if st == -1:
            w = rc(w)
            off = (hi - (p + 3))
        else:
            off = p - lo
        fh.write(f">{c}|{p}|{st}|{off}\n{w}\n")   # window in the gene's orientation; codon at [off, off+3)
hits = D / "out" / "crossmap_blastn.tsv"
if not hits.exists():
    subprocess.run(["blastn", "-query", str(q), "-db", str(D / "out" / "blastdb" / "cngb_other"), "-task", "blastn",
                    "-evalue", "1e-10", "-max_target_seqs", "20", "-num_threads", "4", "-outfmt",
                    "6 qseqid sseqid pident length qstart qend sstart send qseq sseq", "-out", str(hits)], check=True)
best = {}
for line in open(hits):
    qid, sid, pid, ln, qs, qe, ss, se, qa, sa = line.rstrip("\n").split("\t")
    c, p, st, off = qid.split("|")
    off = int(off)
    qs, qe, ln, pid = int(qs), int(qe), int(ln), float(pid)
    if not (qs - 1 <= off and off + 3 <= qe):
        continue
    qi, cod = qs - 2, ""
    for a, b in zip(qa, sa):
        if a != "-":
            qi += 1
        if off <= qi < off + 3 and a != "-":
            cod += b
    key = (c, int(p))
    score = (pid, ln)
    if ln >= 100 and (key not in best or score > best[key][0]):
        best[key] = (score, sid.split("_contig")[0] if "_contig" in sid else sid, cod)
rows = []
for k, per in obs.items():
    alts = Counter()
    for run, o in per.items():
        for cod, n in o.items():
            if cod != "TAG" and n >= 2:
                alts[cod] += 1
    b = best.get(k)
    rows.append((k, per, alts, b))
with open(D / "out" / "crossmap.tsv", "w") as fh:
    fh.write("contig\tpos0\tpolymorphic_alts\tbest_other_pid\tlen\tother_bin\tother_codon\tper_run\n")
    for k, per, alts, b in sorted(rows):
        pr = " | ".join(f"{r[-2:]}:" + ",".join(f"{c}{n}" for c, n in o.most_common()) for r, o in per.items())
        fh.write(f"{k[0]}\t{k[1]}\t{','.join(alts) or '-'}\t{b[0][0] if b else '-'}\t{b[0][1] if b else '-'}\t"
                 f"{b[1] if b else '-'}\t{b[2] if b else '-'}\t{pr}\n")
poly = [r for r in rows if r[2]]
non = [r for r in rows if not r[2]]
for name, rs in (("polymorphic", poly), ("non-polymorphic", non)):
    hi = [r for r in rs if r[3] and r[3][0][0] >= 95]
    match = [r for r in hi if r[3][2] in r[2]]
    print(f"{name} TAG sites: {len(rs)}; another species' genome >= 95% identical over >= 100 nt of the window: {len(hi)}"
          f"; of these, the other genome carries the polymorphic alt codon: {len(match)}")
    print(f"   best identity distribution: {Counter(int(r[3][0][0]) // 5 * 5 if r[3] else 'none' for r in rs).most_common()}")
for r in poly:
    print(f"   {r[0][0]}:{r[0][1]} alts={','.join(r[2])} best_other={r[3][0] if r[3] else '-'} {r[3][1] if r[3] else ''} codon={r[3][2] if r[3] else '-'}")
