"""Scan per-sample assemblies (the 49 Ethiopian chicken caecum IDBA assemblies behind the CASERY01 co-assembly) for:
  1. contigs from the CASERY01 population: minimap2 (asm10) against the MAG, >= 95% identity over >= 1 kb;
  2. whether the MAG's in-frame TAG sites are TAG in that chicken too (codon at the aligned position);
  3. every tRNA with a stop anticodon (ARAGORN via pyaragorn), with whether its contig is CASERY01-linked
     (aligns to the MAG anywhere) and the contig's length/GC.

usage: uv run python src/sample_asm_scan.py <asm.fa.gz> [...]
       -> results/samples/<name>.json
"""
import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pyaragorn
import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import ROOT, fetch, read_fasta, revcomp

MAG = "GCA_949290385.1"
OUT = ROOT / "results" / "samples"
STOP_ANTI = {"cta", "tca", "tta", "CTA", "TCA", "TTA"}


def mag_tag_sites():
    """In-frame TAG positions (contig, 0-based position of T, strand) in MAG genes called with TAG as sense."""
    contigs = dict(read_fasta(fetch(MAG)))
    finder = pyrodigal.GeneFinder()
    finder.train(*[s.encode() for s in contigs.values()], translation_table=15)
    sites = []
    for c, s in contigs.items():
        for g in finder.find_genes(s.encode()):
            nt = g.sequence()
            for i in range(0, len(nt) - 3, 3):
                if nt[i:i + 3] == "TAG":
                    if g.strand == 1:
                        sites.append((c, g.begin - 1 + i, 1))
                    else:
                        sites.append((c, g.end - 1 - i, -1))  # position of the T on the minus strand (forward coords)
    return contigs, sites


def main(paths):
    OUT.mkdir(parents=True, exist_ok=True)
    mag, sites = mag_tag_sites()
    with tempfile.TemporaryDirectory() as td:
        ref = Path(td) / "mag.fa"
        ref.write_text("".join(f">{k}\n{v}\n" for k, v in mag.items()))
        for p in paths:
            name = Path(p).name.split(".")[0]
            if (OUT / f"{name}.json").exists():
                continue
            asm = dict(read_fasta(p))
            paf = Path(td) / "a.paf"
            subprocess.run(["minimap2", "-x", "asm10", "-c", "--cs", "-t", "3", str(ref), p, "-o", str(paf)],
                           check=True, stderr=subprocess.DEVNULL)
            linked, aln = set(), []
            for line in open(paf):
                f = line.split("\t")
                qn, ql, qs, qe, strand, tn, tl, ts, te, nm, blk = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
                if blk >= 1000 and nm / blk >= 0.95:
                    linked.add(qn)
                    aln.append((tn, ts, te, qn, qs, qe, strand))
            # TAG sites covered by an alignment: read the sample's codon at the same place (ungapped approximation,
            # checked by requiring the 10 nt around the site to match the MAG exactly)
            same = diff = 0
            diff_codons = {}
            for c, pos, st in sites:
                for tn, ts, te, qn, qs, qe, strand in aln:
                    if tn == c and ts + 10 <= pos <= te - 13:
                        off = pos - ts
                        q = asm[qn]
                        if strand == "+":
                            qp = qs + off
                            win_m, win_q = mag[c][pos - 10:pos + 13], q[qp - 10:qp + 13]
                        else:
                            qp = qe - off - 1
                            win_m, win_q = mag[c][pos - 10:pos + 13], revcomp(q[qp - 12:qp + 11])
                        if len(win_q) != 23 or win_m[:10] != win_q[:10] or win_m[13:] != win_q[13:]:
                            continue
                        cod_m = win_m[10:13] if st == 1 else revcomp(mag[c][pos - 2:pos + 1])
                        cod_q = win_q[10:13] if st == 1 else revcomp(win_q[8:11])
                        if cod_q == cod_m:
                            same += 1
                        else:
                            diff += 1
                            diff_codons[cod_q] = diff_codons.get(cod_q, 0) + 1
                        break
            # stop-anticodon tRNAs anywhere in the sample assembly
            finder = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
            trnas = []
            for k, s in asm.items():
                if len(s) < 80:
                    continue
                for g in finder.find_rna(s.encode()):
                    if g.anticodon in STOP_ANTI or g.amino_acid in ("Pyl", "SeC", "Stop", "Sup", "???"):
                        gc = (s.count("G") + s.count("C")) / len(s)
                        trnas.append(dict(contig=k, len=len(s), gc=round(gc, 3), begin=g.begin, end=g.end,
                                          strand=g.strand, aa=g.amino_acid, anticodon=g.anticodon,
                                          energy=round(g.energy, 1), linked=k in linked, seq=g.sequence()))
            res = dict(sample=name, contigs=len(asm), bp=sum(len(s) for s in asm.values()), linked_contigs=len(linked),
                       linked_bp=sum(len(asm[k]) for k in linked), mag_aligned_bp=sum(te - ts for tn, ts, te, *_ in aln),
                       tag_sites_checked=same + diff, tag_same=same, tag_diff=diff, tag_diff_codons=diff_codons,
                       stop_trnas=trnas)
            (OUT / f"{name}.json").write_text(json.dumps(res, indent=1))
            print(name, {k: v for k, v in res.items() if k != "stop_trnas"}, "stop tRNAs:", len(trnas),
                  "linked:", sum(t["linked"] for t in trnas), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
