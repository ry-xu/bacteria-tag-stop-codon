"""RF1 (prfA) neighbourhood: call genes (table 15, so TAG readers' genes are not cut), annotate by best Swiss-Prot hit,
and print the genes around RF1/RF2 and around the usual prfA neighbours (lepA, topoisomerase III, bL31, PrmC/HemK).
usage: python src/rf1_locus.py <acc> [<acc> ...]  -> results/locus/<acc>.t15.tsv (gene, contig, begin, end, strand, hit)"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import DB, ROOT, fetch, read_fasta

KEYS = {"RF1_": "RF1", "RF2_": "RF2", "LEPA_": "lepA", "TOP3_": "topIII", "RL31_": "bL31", "PRMC_": "prmC",
        "RL31B_": "bL31"}


def annotate(acc):
    out = ROOT / "results" / "locus" / f"{acc}.t15.tsv"
    if out.exists():
        return [l.rstrip("\n").split("\t") for l in open(out)]
    contigs = dict(read_fasta(fetch(acc)))
    f = pyrodigal.GeneFinder()
    f.train(*[s.encode() for s in contigs.values()], translation_table=15)
    genes = []
    for c, s in contigs.items():
        for g in f.find_genes(s.encode()):
            genes.append((f"g{len(genes)}", c, g.begin, g.end, g.strand, g.translate(translation_table=15)))
    with tempfile.TemporaryDirectory() as td:
        q, o = Path(td) / "q.faa", Path(td) / "o.tsv"
        q.write_text("".join(f">{g[0]}\n{g[5].rstrip('*')}\n" for g in genes))
        subprocess.run(["diamond", "blastp", "-q", q, "-d", DB, "-o", o, "--quiet", "--max-target-seqs", "1",
                        "--evalue", "1e-5", "--threads", "4", "--outfmt", "6", "qseqid", "pident", "sseqid"], check=True)
        hit = {}
        for l in open(o):
            a, p, sid = l.rstrip("\n").split("\t")
            hit.setdefault(a, f"{sid.split('|')[2]}:{p}")
    rows = [[g[0], g[1], str(g[2]), str(g[3]), str(g[4]), hit.get(g[0], "-")] for g in genes]
    out.parent.mkdir(exist_ok=True)
    out.write_text("".join("\t".join(r) + "\n" for r in rows))
    return rows


for acc in sys.argv[1:]:
    rows = annotate(acc)
    print(f"== {acc}")
    idx = [i for i, r in enumerate(rows) if any(r[5].startswith(k) or f"_{k}" in r[5] or r[5].split(":")[0].split("_")[0] + "_" == k
                                               for k in KEYS)]
    shown = set()
    for i in idx:
        for j in range(max(0, i - 2), min(len(rows), i + 3)):
            if j in shown or rows[j][1] != rows[i][1]:
                continue
            shown.add(j)
            r = rows[j]
            tag = KEYS.get(r[5].split(":")[0].split("_")[0] + "_", "")
            print(f"  {r[1]:<22} {r[2]:>8} {r[3]:>8} {'+' if r[4] == '1' else '-'}  {r[5]:<28} {tag}")
        print("  ..")
